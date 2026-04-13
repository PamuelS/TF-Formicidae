SAMPLES = glob_wildcards("OG_used/{samples}.fa")[0]

rule all:
	input:
		expand("../04_trimmed/{samples}_trimmed.faa", samples=SAMPLES)
rule aligned:
	input:
		"OG_used/{samples}.fa"
	output:
		"{samples}_aligned.faa"
	conda:
		"sequence"
	shell:
		"mafft --auto --anysymbol {input} > {output}"
rule trimmed:
	input:
		"{samples}_aligned.faa"
	output:
		trim = "../04_trimmed/{samples}_trimmed.faa",
		html = "../04_trimmed/{samples}"
	conda:
		"sequence"
	shell:
		"bmge -i {input} -of {output.trim} -oh {output.html} -t AA -h 0.5 -g 0.4"
