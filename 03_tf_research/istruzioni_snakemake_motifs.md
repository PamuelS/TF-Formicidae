# Istruzioni
Guida per le varie rule del file di snakemake ottimizzato per lo studio di tutti i motivi associati alle varie specie.

## Creazione Orthogroups_DISCO.tsv
Essendo che nella nostra pipeline di lavoro abbiamo eseguito DISCO per eliminare tutti gli ortogruppi contenenti geni paraloghi, si è optato per la generazione di un file contenent e la proteina associata ad ogni specie che ritroviamo in un ortogruppo specifico (analogo del file `Orthogroups.tsv`)
```bash
awk 'BEGIN {FS = "|"; OFS = "\t"} FNR == 1 {og = FILENAME; gsub(".*/", "", og); gsub(".faa", "", og)} /^>/ {split(substr($0, 2), parts, "|"); sp = parts[1]; gene = parts[2]; data[og, sp] = gene; ogs[og] = 1; spp[sp] = 1} END {header = "Orthogroup"; for (sp in spp) {header = header OFS sp; col_order[++n] = sp} print header; for (og in ogs) {row = og; for (i = 1; i <= n; i++) {sp = col_order[i]; if ((og, sp) in data) {row = row OFS data[og, sp]} else {row = row OFS ""}} print row}}' *.faa > Orthogroups_DISCO.tsv
```

## Snakemake's rules
In totale per eseguire la pipeline di lavoro completa sono state eseguite 12 differenti rule di snakemake.
