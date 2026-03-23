SAMPLES = glob_wildcards("00_genome/{sample}.fna")[0]

rule all:
	input:
		expand("{sample}_busco", sample=SAMPLES)
rule busco:
	input:
		"00_genome/{sample}.fna"
	output:
		directory("{sample}_busco")
	shell:
		"busco -m geno -l /DATASMALL/samuel.pederzini/TF-Formicidae/01_BUSCO/busco_downloads/lineages/hymenoptera_odb10 -c 8 -o {output} -i {input}"
