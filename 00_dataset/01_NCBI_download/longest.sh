rule all:
	input:
		expand("{sample}_longest.faa", sample=SAMPLES)

rule longest_isoform:
	input:
		"{sample}.gff"
	outpu:
		"{sample}_longest.gff"
	shell:
		"agat_sp_keep_longest_isoform.pl --gff {sample}.gff -o {samples}_longest.gff "

rule cds_extraction:
	input:
		gff = "{sample}_longest.gff"
		fasta = "../00_genome/{sample}.fna"
	output:
		"{sample}.faa"
	shell:
		"agat_sp_extract_sequences.pl -g {input.gff} -f ../00_genome/{input.fasta} -t cds -p --cfs --output ../02_Proteome"

