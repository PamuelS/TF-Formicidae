# Allineamento delle sequenze
In questa cartella vengono svolte tutte le operazioni necessarie ad eseguire gli allineamenti delle sequenze dei geni ortologhi single-copy ottenuti dai passaggi precedenti.
Per raggiungere un numero di ortogruppi pari a 200, sono stati utilizzati tutti i `single-copy` ottenuti dall'[analisi di OrthoFinder](../00_Orthofinder_analysis/Statistics_Overall.tsv) e per raggiungere il numero di ortogruppi prestabilito sono stati scelti gli ortogruppi originatisi dopo il lancio delcomando `alternative_split_disco_outgroup.sh`

I mancanti 71 ortogruppi `single-copy` sono stati scelti randomicamente dai risultatti di DISCO
```bash
less every_species_OGs.txt | shuf -n 71 > species_tree.txt
```

> il file `every_species_OGs.txt` indica quanti sono gli ortogruppi che risultanti come single-copy al termine di `alternative_split_disco_outgroup.sh`

successivamente sono stati salvati tutti gli OG utilizzate per le successive analisi in un unico file. Quindi il file `OG_tree.txt` contiene esattamente 200 ortogruppi caratteriuzzati ciascuno da 175 header (ovvero specie)
```bash
ls > ../../02_DISCO_OG/OG_tree.txt
```
> il file lo ritrovi [qui](../02_DISCO_OG/OG_tree.txt)


Per verificare se gli ortogruppi presi randomicamente non si ripetessero con i `single-copy` di OrthoFinder, è stato eseguito il successivo comando:
```bash
ls | sort | uniq -c | less
```

## Analisi
Una volta trovati quelli che sono gli ortogruppi da utilizzare per l'analisi, è stato eseguito uno script snakemake `snakemake_alligned_trimmed.smk` che nel complesso ha eseguito primal'allineamneto delle sequenze (prima rule) e secondariamente lo stesso script ha eseguito anche il triming delle sequenze (seconda rule).
```bash
snakemake -s snakemake_alligned_trimmed.smk --use-conda --cores 5
```
