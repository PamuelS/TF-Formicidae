SAMPLES = glob_wildcards("every_genome/{samples}.fna")[0]
MOTIFS = glob_wildcards("../00_download_motif/{motifs}.jaspar")[0]
BOWTIE_SUFFIX = ["1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"]

rule all:
	input:
		expand("04_bowtie/02_genome_motif_table/{samples}_summary.tsv", samples=SAMPLES),
		expand("04_bowtie/06_full_promoter_motif_table/{samples}_full_summary.tsv", samples=SAMPLES),
		expand("05_aggregate/02_totalscore/totalscore_{motif}.tsv", motif=MOTIFS)

# Creazione delle sequenze promotrici definite come 5k pb upstream e 2k pb downstream
rule create_pep_fasta:
	input:
		genome = "every_genome/{samples}.fna",
		gff = "every_gff/{samples}_longest.gff"
	output:
		fai = temp("02_genome_analysis/00_genome_coordinates/{samples}.fai"),
		size = "02_genome_analysis/01_genome_size/{samples}.size",
		tss = "02_genome_analysis/02_tss/{samples}_tss.bed",
		bed = "02_genome_analysis/03_promoter_bed/{samples}.bed",
		fasta = "02_genome_analysis/04_promoter_5k_2k_fasta/{samples}_5k_2k.fasta"
	conda:
		"bedtool"
	shell:
		"""
		samtools faidx --fai-idx {output.fai} {input.genome}
		cut -f1,2 {output.fai} > {output.size}
		awk -v FS="\t" '$3=="mRNA" {{split($9, a, "ID="); split(a[2], b, ";"); geneid=b[1]; if($7=="+") print $1"\t"$4-1"\t"$4"\t"geneid"\t.\t"$7; else print $1"\t"$5-1"\t"$5"\t"geneid"\t.\t"$7}}' {input.gff} > {output.tss}
		bedtools slop -i {output.tss} -g {output.size} -l 5000 -r 2000 -s > {output.bed}
		bedtools getfasta -fi {input.genome} -bed {output.bed} -s -name -fo {output.fasta}
		"""

# Indicizzazione per i promotori/peptidi
rule build_promoter_indices:
	input:
		pep_fastas = "02_genome_analysis/04_promoter_5k_2k_fasta/{samples}_5k_2k.fasta"
	output:
		indices = expand("04_bowtie/03_promoter_indices/{{samples}}.{suffix}", suffix=BOWTIE_SUFFIX)
	conda:
		"jaspar"
	params:
		prefix = "04_bowtie/03_promoter_indices/{samples}"
	resources:
		mem=16000,
		time=120,
	shell:
		"bowtie-build -f {input.pep_fastas} {params.prefix}"

# Indicizzazione di tutti i genomi completi
rule genomes_indices:
	input:
		genome = "every_genome/{samples}.fna"
	output:
		indices = expand("04_bowtie/00_genome_indices/{{samples}}.{suffix}", suffix=BOWTIE_SUFFIX)
	conda:
		"jaspar"
	params:
		basename = "04_bowtie/00_genome_indices/{samples}"
	shell:
		"bowtie-build -f {input.genome} {params.basename}"

# Creazione del background
rule background_comp:
	input:
		"every_genome/{samples}.fna"
	output:
		"01_background_genomes/{samples}_bg.tsv"
	shell:
		"seq_extract_bcomp -i 0 -c {input} | "
		"awk -F, '{{print \"A\\t\" $1 \"\\nC\\t\" $2 \"\\nG\\t\" $3 \"\\nT\\t\" $4}}' > {output}"

# Conversione matrici PFM in PWM
rule pwm_conv:
	input:
		motifs = "00_download_motif/{motifs}.jaspar",
		bg = "01_background_genomes/{samples}_bg.tsv"
	output:
		ill = "03_conversion/00_pwm_convert/{samples}/{samples}_{motifs}.ill",
		score = "03_conversion/01_pwm_score/{samples}/{samples}_{motifs}.score",
		mba = "03_conversion/02_pwm_mba/{samples}/{samples}_{motifs}.mba",
		tags = "03_conversion/03_tags/{samples}/{samples}_{motifs}.tags"
	threads:
		1
	shell:
		"""
		pwm_convert {input.motifs} -f=jaspar -b={input.bg} > {output.ill}
		BG_VALUES=$(awk '{{print $2}}' {input.bg} | grep -E '^[0-9.]+$' | tr '\\n' ',' | sed 's/,$//')
		matrix_prob -e 0.00001 --bg $BG_VALUES {output.ill} > {output.score}
		SCORE=$(grep 'SCORE :' {output.score} | awk '{{print $3}}')
		mba -c $SCORE {output.ill} > {output.mba}
		awk '{{print ">"$2"\\n"$1}}' {output.mba} > {output.tags}
		"""

# Mappatura dei tags ai genomi completi — un SAM per motivo
rule genome_map_tags:
	input:
		tags = "03_conversion/03_tags/{samples}/{samples}_{motifs}.tags",
		gen_ind = expand("04_bowtie/00_genome_indices/{{samples}}.{suffix}", suffix=BOWTIE_SUFFIX)
	output:
		mapped_tags = "04_bowtie/01_genome_map_tags/{samples}_{motifs}.sam"
	conda:
		"jaspar"
	params:
		prefix = "04_bowtie/00_genome_indices/{samples}"
	threads:
		1
	resources:
		mem=16000,
		time=120,
	shell:
		"bowtie -n 0 -a {params.prefix} -f {input.tags} > {output.mapped_tags}"

# Mappatura dei tags sui promotori/peptidi — un SAM per motivo
rule promoter_map_tags:
	input:
		tags = "03_conversion/03_tags/{samples}/{samples}_{motifs}.tags",
		indices = expand("04_bowtie/03_promoter_indices/{{samples}}.{suffix}", suffix=BOWTIE_SUFFIX)
	output:
		promoter_tags = "04_bowtie/04_promoter_map_tags/{samples}_{motifs}.sam"
	conda:
		"jaspar"
	params:
		prefix = "04_bowtie/03_promoter_indices/{samples}"
	threads: 1
	resources:
		mem=16000,
		time=120,
	shell:
		"bowtie -n 0 -a {params.prefix} -f {input.tags} > {output.promoter_tags}"

# Tabelle dei motivi per il genoma completo
rule genome_motif_tables:
	input:
		genome_tags = expand("04_bowtie/01_genome_map_tags/{{samples}}_{motifs}.sam", motifs=MOTIFS)
	output:
		genome_species_table = "04_bowtie/02_genome_motif_table/{samples}_summary.tsv"
	threads: 1
	resources:
		mem=2000,
		time=60,
	run:
		from pathlib import Path
		import pandas as pd
		import sys

		out = []
		for motif, file in zip(MOTIFS, input.genome_tags):
			data = pd.read_csv(
				file, sep=r'\s+', header=None, # read in a space separated file without a header
				usecols=[0, 2], #only take the score and scaffold
				names=['score', 'scaffold']) # name them
				#data has a single file's data that will be added to output
			out.append(
				#group the data by scaffold name, select score, aggregate the mean and count of score, reset index to make the scaffold a column, then save the motif name for later
				data.groupby('scaffold')['score'].agg(['mean', 'count']).reset_index().assign(motif=motif)
			)
		# combine all the data tables into one and save as a tab separated file
		pd.concat(out).to_csv(output.genome_species_table, sep='\t', index=False)

# Tabelle dei motivi per i promotori/peptidi
rule species_motif_tables:
	input:
		promoter_tags = expand("04_bowtie/04_promoter_map_tags/{{samples}}_{motifs}.sam", motifs=MOTIFS)
	output:
		species_table = "04_bowtie/05_promoter_motif_table/{samples}_summary.tsv"
	threads: 1
	resources:
		mem=2000,
		time=60,
	run:
		from pathlib import Path
		import pandas as pd
		import sys

		out = []
		for motif, file in zip(MOTIFS, input.promoter_tags):
			data = pd.read_csv(
				file, sep=r'\s+', header=None, # read in a space separated file without a header
				usecols=[0, 2], #only take the score and peptide
				names=['score', 'peptide']) # name them
				#data has a single file's data that will be added to output
			out.append(
				#group the data by peptide name, select score, aggregate the mean and count of score, reset index to make the peptide a column, then save the filename for later
				data.groupby('peptide')['score'].agg(['mean', 'count']).reset_index().assign(motif=motif)
			)
		# combine all the data tables into one and save as a tab separated file
		pd.concat(out).to_csv(output.species_table, sep='\t', index=False)

# Creazione di una species_motif_table che contenga ogni singola proteina estratta dal gff (isoforma più lunga)
rule full_species_motif_tables:
	input:
		species_table = "04_bowtie/05_promoter_motif_table/{samples}_summary.tsv",
		abbreviative = "../00_dataset/00_GAGA_download/GAGA_vs_personal_ID.tsv",
	output:
		full_name = "04_bowtie/06_full_promoter_motif_table/{samples}_full_summary.tsv"
	shell:
		"""
		MAP_PATH="../ncbi_header_modified/{wildcards.samples}_header_summary_mapping.tsv"

		if [ -f "$MAP_PATH" ]; then 
			sed -E "s/::.[^\t]+//" {input.species_table} | \
			awk 'BEGIN {{ OFS="\t" }}
				{{ sub(/\\r$/, "") }} 
				NR==FNR {{ map[$2]=$1; next }} 
				$1 in map {{ $1=map[$1] }} 
				1' "$MAP_PATH" - > {output.full_name}

		else 
			abb=$(awk -v fn="{wildcards.samples}" '$3 == fn {{print $2}}' {input.abbreviative})
			sed -E "s/${{abb}}_?//; s/^_//; s/__/_/g; s/::.[^\t]+//" {input.species_table} > {output.full_name}
		fi
		"""

# Creazione di una species_motif_table che contenga unicamente le proteine esistenti dentro il file Orthogroups_DISCO.tsv
rule disco_species_motif_table:
	input:
		species_table = "04_bowtie/05_promoter_motif_table/{samples}_summary.tsv",
		abbreviative = "/DATASMALL/samuel.pederzini/TF-Formicidae/00_dataset/00_GAGA_download/GAGA_vs_personal_ID.tsv",
		orthogroups = "../Orthogroups_DISCO.tsv"
	output:
		disco_name = "04_bowtie/07_disco_promoter_motif_table/{samples}_disco_summary.tsv"
	shell:
		"""
		MAP_PATH="../ncbi_header_modified/{wildcards.samples}_header_summary_mapping.tsv"

		if [ -f "$MAP_PATH" ]; then
			sed -E "s/::.[^\t]+//" {input.species_table} | \
			awk 'BEGIN {{ OFS="\t" }} {{ sub(/\\r$/, "") }} NR==FNR {{ map[$2]=$1; next }} $1 in map {{ $1=map[$1] }} 1' "$MAP_PATH" - | \
			awk -F'\t' -v species="{wildcards.samples}" '
				BEGIN {{ OFS="\t" }}
				{{ gsub(/\\r$/, "") }}
				NR==FNR {{
					if (FNR == 1) {{
						for (i=1; i<=NF; i++) {{
							if ($i == species) {{ col = i; break; }}
						}}
						next;
					}}
					if (col > 0 && $col != "" && $col != "*") {{
						n = split($col, a, /,[ ]*/);
						for (j=1; j<=n; j++) if (a[j] != "") valid[a[j]] = 1;
					}}
					next;
				}}
				FNR == 1 {{ print; next }}
				$1 in valid
			' {input.orthogroups} - > {output.disco_name}

		else
			abb=$(awk -v fn="{wildcards.samples}" '$3 == fn {{print $2}}' {input.abbreviative})
			sed -E "s/${{abb}}_?//; s/^_//; s/__/_/g; s/::.[^\t]+//" {input.species_table} | \
			awk -F'\t' -v species="{wildcards.samples}" '
				BEGIN {{ OFS="\t" }}
				{{ gsub(/\\r$/, "") }}
				NR==FNR {{
					if (FNR == 1) {{
						for (i=1; i<=NF; i++) {{
							if ($i == species) {{ col = i; break; }}
						}}
						next;
					}}
					if (col > 0 && $col != "" && $col != "*") {{
						n = split($col, a, /,[ ]*/);
						for (j=1; j<=n; j++) if (a[j] != "") valid[a[j]] = 1;
					}}
					next;
				}}
				FNR == 1 {{ print; next }}
				$1 in valid
			' {input.orthogroups} - > {output.disco_name}
		fi
		"""
# Aggregazione finale delle tabelle
rule aggregate_tables:
	input:
		species_tables = expand("04_bowtie/07_disco_promoter_motif_table/{samples}_disco_summary.tsv", samples=SAMPLES),
		orthogroups = "../Orthogroups_DISCO.tsv"
	output:
		score_tables = expand("05_aggregate/00_score/score_{motif}.tsv", motif=MOTIFS),
		count_tables = expand("05_aggregate/01_count/count_{motif}.tsv", motif=MOTIFS),
		totalscore_tables = expand("05_aggregate/02_totalscore/totalscore_{motif}.tsv", motif=MOTIFS)
	threads: 1
	resources:
		mem=8000,
		time=120,
	run:
		# Insieme di pacchetti scaricati per eseguire le successive operazioni
		from pathlib import Path
		import sys
		import pandas as pd 
		from collections import defaultdict

		# Creazione del dizionario di traduzione
		# Ricerca il nome del peptide all'interno del file Orthogroups_DISCO.tsv
		# La funzione legge l'header dell'ortogruppo e trova l'indice della colonna corrispondente a quel nome
		# Per ogni riga del file estrae il nome dell'OG e la lunga stringa di nomi di proteine separata da virgole 
		# Fa in modo che ogni nome di peptidfe punti direttamente al suo OG di riferimento
		def build_orthogroup_dict(file, species):
			result = defaultdict(lambda: 'NA')
			header = file.readline()
			spp = header.strip().split('\t')
			index = -1
			for i, sp in enumerate(spp):
				if sp.startswith(species):
					index = i
					break
			else:
				raise ValueError(f"Species {species} not found in header")
			for line in file:
				columns = line.strip('\n').split('\t')
				orthogroup, peptides = columns[0], columns[index]
				peptides = peptides.split(',')
				result.update({peptide.strip(): orthogroup for peptide in peptides})
			return result

		scores = defaultdict(list)
		counts = defaultdict(list)
		totalscores = defaultdict(list)

		# Lo script scorre ogni nome di specie e il relativo file "_summary.tsv"
		for sample, file in zip(SAMPLES, input.species_table):
			print(sample)
			# Per ciascuna specie crea il dizionario proteina-OG
			db = build_orthogroup_dict(open(input.orthogroups), sample)
			# rimozione keys vuote (ortogruppo no proteina nella specie)
			del db['']
			# rimozione keys con values NA (ad esempio header)
			db = {k: v for k, v in db.items() if not pd.isna(v) and v != "NA"}
			# Lettura DISCO_summary
			data = pd.read_csv(file, sep=r'\s+')

			
			# Mappa ogni singolo peptide (colonna di data) con il OG di riferimento (presente in db)
			data['OG'] = data.peptide.map(db)
			# Calcolo complessivo dei punteggi (media * conteggio) per ciascun OG
			data['total_score'] = data['mean'] * data['count']
			for motif in MOTIFS:

				# Qui lo script inizia a girare per ogni singolo motivo
				# Isola le righe di dati che corrispondono solo al motivo target 
				# Comprime tutti i dati ottenuti (mean && count) dentro il singolo nome di OG
				# Calcola il complessivo risultato di per ogni OG (total_count)
				sub_data = data[data.motif == motif]
				# filtra la tabella
				sub_data = sub_data[["OG", "mean", "count", "total_score"]]
				# imposta OG come rowname (index)
				sub_data = sub_data.set_index("OG")
				# rimuovi la riga che adesso ha come nome il nome della colonna ed è quindi vuota
				sub_data.index.name = None

				sub_data = pd.concat((
					sub_data, 
					pd.DataFrame(
						index=sorted(set(db.values()).difference(sub_data.index)),
						columns=sub_data.columns, 
						data=0
					)
				)).sort_index()
				scores[motif].append(sub_data['mean'].rename(sample))
				counts[motif].append(sub_data['count'].rename(sample))
				totalscores[motif].append(sub_data['total_score'].rename(sample))

		# Una volta prodotti i dati necessari, vengono creati tutti i tre file necessari per la lettura dei risultati
		for motif, score, count, totalscore in zip(MOTIFS, output.score_tables, output.count_tables, output.totalscore_tables):
			pd.concat(scores[motif], axis=1, ignore_index=False).reset_index().to_csv(score, sep='\t', index=False, na_rep='NA')
			pd.concat(counts[motif], axis=1, ignore_index=False).reset_index().to_csv(count, sep='\t', index=False, na_rep='NA')
			pd.concat(totalscores[motif], axis=1, ignore_index=False).reset_index().to_csv(totalscore, sep='\t', index=False, na_rep='NA')
