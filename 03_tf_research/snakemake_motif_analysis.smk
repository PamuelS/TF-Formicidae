SAMPLES = glob_wildcards("every_genome/{samples}.fna")[0]
MOTIFS = glob_wildcards("00_download_motif/{motifs}.jaspar")[0]
BOWTIE_SUFFIX = ["1.ebwt", "2.ebwt", "3.ebwt", "4.ebwt", "rev.1.ebwt", "rev.2.ebwt"]

rule all:
	input:
		expand("04_bowtie/01_genome_map_tags/{samples}.sam", samples=SAMPLES),
		expand("05_aggregate/02_totalcount/totalcount_{motif}.tsv", samples=SAMPLES, motif=MOTIFS)

# Creazione dei file che contengono solamente due colonne (abbreviativo_specie   nome_proteina)
rule create_species_pep_pairs:
	input:
		gff = "every_gff/{samples}_representative.gff3"
	output:
		pairs = "02_genome_analysis/00_species_pep_pairs/{samples}_pep_pairs.tsv"
	threads: 1
	resources:
		mem=2000,
		time=60,
	run:
		import pandas as pd

		df = pd.read_csv(input.gff, comment='#', sep='\t', header=None, low_memory=False)
		mRNAs = df[df[2] == 'mRNA'].copy()

		final_data = []
		if not mRNAs.empty:
			mRNAs['length'] = mRNAs[4] - mRNAs[3]
			mRNAs['peptide_id'] = mRNAs[8].str.extract(r'ID=([^;]+)')
			mRNAs['parent_gene'] = mRNAs[8].str.extract(r'Parent=([^;]+)')
			longest = mRNAs.sort_values('length', ascending=False).drop_duplicates(subset=['parent_gene'])
			for pep_id in longest['peptide_id'].dropna():
				final_data.append([wildcards.samples, pep_id])

		pd.DataFrame(final_data).to_csv(output.pairs, sep='\t', header=False, index=False)


# Creazione dei file BED per ciascuna delle specie
rule create_pep_beds:
	input:
		isoforms = "every_gff/{samples}_representative.gff3",
		species_peptide_pairs = "02_genome_analysis/00_species_pep_pairs/{samples}_pep_pairs.tsv"
	output:
		"02_genome_analysis/01_bed_files/{samples}.bed"
	threads: 1
	resources:
		mem=2000,
		time=60,
	run:
		import pandas as pd
		species_peptides = pd.read_csv(input.species_peptide_pairs, sep='\t', names=['species', 'peptide'])
		species_peptides = species_peptides[species_peptides['species'] == wildcards.samples]
		isoforms = pd.read_csv(input.isoforms, comment='#', sep='\t', header=None)
		isoforms = isoforms[isoforms[2] == 'mRNA']
		isoforms['peptide'] = isoforms[8].str.extract(r'ID=([^;]+)')[0]
		isoforms['dot'] = '.'
		isoforms = isoforms[isoforms['peptide'].isin(species_peptides['peptide'])]
		isoforms[[0, 3, 4, 'peptide', 'dot', 6]].to_csv(output[0], sep='\t', header=None, index=None)

# Creazione dei file FASTA per ciascuna delle specie
rule create_pep_fastas:
	input:
		bed = "02_genome_analysis/01_bed_files/{samples}.bed",
		genome = "every_genome/{samples}.fna"
	conda:
		"bedtool"
	output:
		pep_fastas = "02_genome_analysis/02_pep_fasta/{samples}.fasta"
	threads: 1
	resources:
		mem=2000,
		time=60,
	run:
		import pybedtools
		from Bio import SeqIO
		from Bio.SeqRecord import SeqRecord
		from Bio.Seq import Seq
		bed = pybedtools.BedTool(input.bed)
		genome = SeqIO.to_dict(SeqIO.parse(input.genome, "fasta"))
		records = []
		for region in bed:
			chrom = region.chrom
			start = int(region.start)
			end = int(region.end)
			name = region.name
			seq = genome[chrom].seq[start:end]
			record = SeqRecord(Seq(str(seq)), id=name, description=chrom+":"+str(region.start)+"-"+str(region.end))
			records.append(record)
		SeqIO.write(records, output.pep_fastas, "fasta")

# Indicizzazione per i promotori/peptidi
# ricordarsi che l'input della rule sottostante è stata riprodotta basandosi sulla pipeline di snakemake dei TF delle api (quindi non è stato utilizzata la cartella 02_genome_analysis/06_promoters_fasta)
rule build_promoter_indices:
	input:
		pep_fastas = "02_genome_analysis/02_pep_fasta/{samples}.fasta"
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
		"every_genome/{genomes}.fna"
	output:
		"01_background_genomes/{genomes}_bg.tsv"
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

# Mappatura dei tags ai genomi completi
rule genome_map_tags:
	input:
		taglist = expand("03_conversion/03_tags/{samples}/{samples}_{motifs}.tags", samples=SAMPLES, motifs=MOTIFS),
		gen_ind = expand("04_bowtie/00_genome_indices/{{samples}}.{suffix}", suffix=BOWTIE_SUFFIX)
	output:
		mapped_tags = "04_bowtie/01_genome_map_tags/{samples}.sam"
	conda:
		"jaspar"
	params:
		prefix = "04_bowtie/03_promoter_indices/{samples}",
		tag_string = lambda wildcards, input: ",".join(input.taglist)
	threads:
		10
	shell:
		"bowtie -p {threads} -n 0 -a --sam {params.prefix} -f {params.tag_string} > {output.mapped_tags}"

# Mappatura dei tags sui promotori/peptidi
rule promoter_map_tags:
	input:
		taglist = expand("03_conversion/03_tags/{samples}/{samples}_{motifs}.tags", samples=SAMPLES, motifs=MOTIFS),
		indices = expand("04_bowtie/03_promoter_indices/{{samples}}.{suffix}", suffix=BOWTIE_SUFFIX)
	output:
		promoter_tags = "04_bowtie/04_promoter_map_tags/{samples}.sam"
	conda:
		"jaspar"
	params:
		prefix = "04_bowtie/03_promoter_indices/{samples}",
		tag_string = lambda wildcards, input: ",".join(input.taglist)
	threads: 10
	resources:
		mem=16000,
		time=120,
	shell:
		"bowtie -p {threads} -n 0 -a {params.prefix} -f {params.tag_string} > {output.promoter_tags}"

# Tabelle dei motivi per il genoma completo
rule genome_motif_tables:
	input:
		mapped_tags = "04_bowtie/01_genome_map_tags/{samples}.sam"
	output:
		genome_species_table = "04_bowtie/02_genome_motif_table/{samples}_summary.tsv"
	threads:
		20
	run:    
		import pandas as pd
		import sys
		data = pd.read_csv(input.mapped_tags, sep='\t', comment='@', header=None, usecols=[0, 2, 4], names=['motif', 'scaffold', 'score'])
		summary = (data.groupby(['motif', 'scaffold'])['score'].agg(['mean', 'count']).reset_index())
		summary.to_csv(output.genome_species_table, sep='\t', index=False)

# Tabelle dei motivi per i promotori/peptidi
rule species_motif_tables:
	input:
		promoter_tags = "04_bowtie/04_promoter_map_tags/{samples}.sam"
	output:
		species_table = "04_bowtie/05_promoter_motif_table/{samples}_summary.tsv"
	threads: 20
	resources:
		mem=2000,
		time=60,
	run:
		import pandas as pd
		import sys
		data = pd.read_csv(input.promoter_tags, sep='\t', comment='@', header=None, usecols=[0, 2, 4], names=['motif', 'peptide', 'score'])
		summary = (data.groupby(['motif', 'peptide'])['score'].agg(['mean', 'count']).reset_index())
		summary.to_csv(output.species_table, sep='\t', index=False)

# Aggregazione finale delle tabelle
rule aggregate_tables:
	input:
		species_tables = expand("04_bowtie/05_promoter_motif_table/{samples}_summary.tsv", samples=SAMPLES),
		orthogroups = "Orthogroups_DISCO.tsv",
		map = "/DATASMALL/samuel.pederzini/TF-Formicidae/00_dataset/00_GAGA_download/GAGA_vs_personal_ID.tsv"
	output:
		score_tables = expand("05_aggregate/00_score/score_{motif}.tsv", samples=SAMPLES, motif=MOTIFS),
		count_tables = expand("05_aggregate/01_count/count_{motif}.tsv", samples=SAMPLES, motif=MOTIFS),
		totalcount_tables = expand("05_aggregate/02_totalcount/totalcount_{motif}.tsv", samples=SAMPLES, motif=MOTIFS)
	threads: 1
	resources:
		mem=8000,
		time=120,
	run:
		from pathlib import Path
		import os
		import sys
		import pandas as pd
		from collections import defaultdict
		import re

		# Caricamento mapping sample -> species_short e species_full
		sample_to_species_short = {}
		sample_to_species_full = {}
		with open(input.map) as f:
			for line in f:
				parts = line.strip().split('\t')
				sample_id, species_short, species_full = parts[0], parts[1], parts[2]
				sample_to_species_short[species_full] = species_short  # Acafer -> Afer
				sample_to_species_full[species_full] = species_full    # Acafer -> Acafer

		def normalize_peptide(name, species_short):
			"""Rimuove species_short dal nome del peptide ovunque si trovi"""
			normalized = re.sub(rf'^{re.escape(species_short)}_', '', name)
			normalized = re.sub(rf'_?{re.escape(species_short)}_?', '_', normalized)
			normalized = re.sub(r'_+', '_', normalized).strip('_')
			return normalized

		def build_orthogroup_dict(file, species_full):
			"""Costruisce {gene_id -> ortogruppo} leggendo Orthogroups_DISCO.tsv"""
			result = defaultdict(lambda: 'NA')
			header = file.readline()
			spp = header.strip().split('\t')
			index = -1
			for i, sp in enumerate(spp):
				if sp == species_full:
					index = i
					break
			else:
				raise ValueError(f"Species {species_full} not found in header")
			for line in file:
				columns = line.strip('\n').split('\t')
				orthogroup = columns[0]
				if index >= len(columns) or columns[index] == '':
					continue
				gene = columns[index].strip()
				result[gene] = orthogroup
			return result

		scores = defaultdict(list)
		counts = defaultdict(list)
		totalcounts = defaultdict(list)

		for sample, file in zip(SAMPLES, input.species_tables):
			print(sample)
			species_short = sample_to_species_short[sample]  # Afer
			species_full = sample_to_species_full[sample]    # Acafer

			db = build_orthogroup_dict(open(input.orthogroups), species_full)
			if '' in db:
				del db['']

			data = pd.read_csv(file, sep='\t')
			data['peptide'] = data['peptide'].map(lambda x: normalize_peptide(x, species_short))
			data['OG'] = data['peptide'].map(db)
			data['total_score'] = data['mean'] * data['count']

			# Diagnostica
			n_na = (data['OG'] == 'NA').sum()
			n_total = len(data)
			print(f"{sample}: {n_na}/{n_total} peptidi non mappati ({100*n_na/n_total:.1f}%)")
			if n_na > 0:
				print("Esempi non mappati:")
				print(data[data['OG'] == 'NA']['peptide'].head(10).tolist())

			for motif in MOTIFS:
				sub_data = data[data.motif == motif]
				summed_data = sub_data.groupby('OG').agg(
					total_score=('total_score', 'sum'),
					total_count=('count', 'sum'),
					mean_count=('count', 'mean')
				)
				summed_data['mean_score'] = summed_data['total_score'] / summed_data['total_count']
				summed_data = pd.concat((summed_data, pd.DataFrame(
					index=sorted(set(db.values()).difference(summed_data.index)),
					columns=summed_data.columns, data=0))).sort_index()

				scores[motif].append(summed_data['mean_score'].rename(sample))
				counts[motif].append(summed_data['mean_count'].rename(sample))
				totalcounts[motif].append(summed_data['total_count'].rename(sample))

		for motif, score_file, count_file, totalcount_file in zip(MOTIFS, output.score_tables, output.count_tables, output.totalcount_tables):
			pd.concat(scores[motif], axis=1, ignore_index=False).reset_index().to_csv(score_file, sep='\t', index=False, na_rep='NA')
			pd.concat(counts[motif], axis=1, ignore_index=False).reset_index().to_csv(count_file, sep='\t', index=False, na_rep='NA')
			pd.concat(totalcounts[motif], axis=1, ignore_index=False).reset_index().to_csv(totalcount_file, sep='\t', index=False, na_rep='NA')
