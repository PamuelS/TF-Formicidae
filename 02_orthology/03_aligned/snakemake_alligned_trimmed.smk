# script used for the allignment and the trimming for every OGs present after the DISCO launch
# for the 200 OGs one used for the creation of the ML tree, simply change the input and use the folder present in 03_alligned 

SAMPLES = glob_wildcards("../02_DISCO_OG/{samples}.faa")[0]

rule all:
        input:
                expand("../04_trimmed/{samples}_trimmed.faa", samples=SAMPLES)
rule aligned:
        input:
                "../02_DISCO_OG/{samples}.faa"
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
                html = "../04_trimmed/{samples}.html"
        conda:
                "sequence"
        shell:
                "bmge -i {input} -of {output.trim} -oh {output.html} -t AA -h 0.5 -g 0.4"

