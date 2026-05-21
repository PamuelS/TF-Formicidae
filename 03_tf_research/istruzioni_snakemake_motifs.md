# Istruzioni
Guida per le varie rule del file di snakemake ottimizzato per lo studio di tutti i motivi associati alle varie specie.

## Creazione Orthogroups_DISCO.tsv
Essendo che nella nostra pipeline di lavoro abbiamo eseguito DISCO per eliminare tutti gli ortogruppi contenenti geni paraloghi, si è optato per la generazione di un file contenent e la proteina associata ad ogni specie che ritroviamo in un ortogruppo specifico (analogo del file `Orthogroups.tsv`)
```bash
awk 'BEGIN {FS = "|"; OFS = "\t"} FNR == 1 {og = FILENAME; gsub(".*/", "", og); gsub(".faa", "", og)} /^>/ {split(substr($0, 2), parts, "|"); sp = parts[1]; gene = parts[2]; data[og, sp] = gene; ogs[og] = 1; spp[sp] = 1} END {header = "Orthogroup"; for (sp in spp) {header = header OFS sp; col_order[++n] = sp} print header; for (og in ogs) {row = og; for (i = 1; i <= n; i++) {sp = col_order[i]; if ((og, sp) in data) {row = row OFS data[og, sp]} else {row = row OFS ""}} print row}}' *.faa > Orthogroups_DISCO.tsv
```

## Lancio di Snakemake
Lo [script di snakemake](./snakemake_motif_study.smk) è stato lanciato second oi lseguente comando

```bash
snakemake -s snakemake_motif_study.smk --use-conda --cores 35 --rerun-incomplete --touch
```

## Snakemake's rules
In totale per eseguire la pipeline di lavoro completa sono state eseguite 12 differenti rule di snakemake, riassunte nelle seguenti informazioni:

***Rule create_pep_fasta***:
Questa rule rappresenta la parte più cruciale per approcciare lo studio dei motivi all'interno delle sequenze promotrici. In questa sezione di codice viene affrontata l'estrazione delle sequenze promotrici mediante alcuni passaggi. tramite l'uso di strumenti bioinformatici (quali samtools e bedtools) vengono estratte le sequenze promotrici, a partire dalle informazioni contenute dentro i file `_longest.gff`, che sono state caratterizzate da sequenze di 5000 pb a monte e 2000 pb a valle delle sequenze TSS (Transcription Starting Site). In seguito all'identificazione di tali porzioni, si associano le sequenze nucleotidiche estratte in questo caso dai file `.fna`.

***Rule build_promoter_indices***:
Questa sezione dello snakemake file consente di andare a creare un indice (suddiviso in sei differenti file) di tutti i promotori che sono stati identificati. Questo procedimento di indicizzazione è stato eseguito con il programma Bowtie v1 (non la v2). 

***Rule genomes_indices***:
Lo stesso procedimento adottato per l’indicizzazione dei promotori, è stato adoperato persino per tutti i genomi. Quindi analogamente al precedente è stato impiegato Bowtie v1.

***Rule background_comp***:
In questa sezione è stato utilizzato il programma seq_extract_bcomp (proveniente da PWMScan) per creare un background di riferimento per ogni specie. Ovvero, viene creato un file contenente l’abbondanza delle varie basi azotate all’interno del genoma espressi in percentuale.

***Rule pwm_convert***:
Questa rule rappresenta una delle più importanti rule adoperate dentro lo script snakemake. La sua funzionalità risiede proprio nella ricostruzioni e adattamento delle sequenze dei motivi (ottenute dal database online JASPAR 2026 CORE) seguendo le caratteristiche del genoma delle specie, mediante l’utilizzo della composizione background delle medesime. L’operazione sostanziale eseguita in questa rule è la conversione della matrice di frequenza in una matrice probabilistica perfettamente adattata alle caratteristiche del genoma analizzato. (Vengono definiti come tags i motivi riadattati alle specie del dataset)

***Rule genome_map_tags***:
Una volta ottenuti tutti e 300 i possibili motivi individualizzati ad hoc per ciascuna specie, si prosegue con la fase di mappatura dei medesimi motivi sul genoma di riferimento. Procedimento eseguito con l’utilizzo del genoma indicizzato e dei tag ovviamente.

***Rule promoter_map_tags***:
Analogamente I tags vengono mappati e collocati sulle sequenze promotrici mediante gli indici dei promotori ottenuti precedentemente. (Procedimento concettualmente identico alla rule precedente)

***Rule genome_motif_tables***:
Viene costruita una tabella a partire dalle informazioni estratte dalla mappatura dei motivi sul genoma e vengono calcolate alcune sattistiche che indicano la compatibilità del motivo con la sequenza associata. La "conta" esprime il quantitativo di volte che quello specifico motivo viene identificato dentro lo scaffold indicato, mentre la "media" consente di capire la robustezza media di quell'appaiamento.

***Rule species_motif_tables***:
In questa rule vengono eseguite le medesime operazione fatte nella precedente rule, con la singola differenza che qui si valuta l'appaiamento del motivo sulle sequenze promotrici e non sul complressivo genoma.

***Rule full_species_motif_tables***:
Quello che viene eseguito in questa rule è la costruzione di una tabella che possegga le nomenclature corrette delle proteine che vengono ripulite del suffisso che le era stato appaiato in precedenza (nella rule create_pep_fasta). Sulla base della specie (se GAGA oppure NCBI) è stato anche intrapreso un tipo di modifica del nome della proteina differenziato:
- per specie GAGA ---> rimozione dell'abbreviativo GAGA in qualunque posizione
- per specie NCBI ---> trasformazione completa del nome basata su un file di riferimento (file contenente nella prima colonna il nome della proteina presente dentro `Orthogroups_DISCO.tsv` e nella seconda colonna il nome della proteina estratto erronemanete da `_longest.gff`)

***Rule disco_species_motif_tables***:
Questa rule esegue sostanzialmente le medesime operazione osservabili nella rule precedente, con la aggiunta della rimozione di ogni singola proteina che non fosse presente e caratterizzata all'interno del file `Orthogroups_DISCO.tsv`.

***Rule aggregate_tables***:
Essendo che vengono utilizzati tutti i file .gff3 provenienti dalle annotazioni GAGA, ritroviamo all'interno una sintassi delle proteine non ancora standardizzata (ovvero con il nome dell'abbreviativo GAGA dentro al nome). Per questo motivo viene eliminato il nome dell'abbreivativo lasciando il semplice nome del peptide, per consentire l'appaiamento corretto con il nome del peptide ritovabile nel Ortogruppo di riferimento (nome già standardizzato precedentemente).
Il risultato finale consiste nella generazione di una singola tabella per ogni motivo associato ad una singola specie, che consente di verificare quali sono gli ortogruppi associati a quello specifico motivo.

<p align="center">
  <img src="./snakemake_motif_study.png" alt="Descrizione">
</p>
