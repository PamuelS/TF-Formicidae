SAMPLES, = glob_wildcards("{sample}.gff")

rule all:
	input:
		expand("{sample}.faa", sample=SAMPLES)

rule longest_isoform:
	input:
		"{sample}.gff"
	output:
		"{sample}_longest.gff"
	conda:
		"sequence"
	shell:
		"agat_sp_keep_longest_isoform.pl --gff {input} -o {output}agat_sp_keep_longest_isoform.pl --gff {input} -o {output} "

rule cds_extraction:
	input:
		gff = "{sample}_longest.gff",
		fasta = "../00_genome/{sample}.fna"
	conda:
                "sequence"
	output:
		"{sample}.faa"
	shell:
		"agat_sp_extract_sequences.pl -g {input.gff} -f ../00_genome/{input.fasta} -t cds -p --cfs --output {output}"

