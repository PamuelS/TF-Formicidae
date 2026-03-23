# Specie NCBI
Per le specie che non erano già presenti all'interno del dataset di GAGA, si è eseguito un approccio differente mediante l'utilizzo di uno [script](./download_dataset.sh) che ha consentito di scaricare sia i genomi (.fna) che le corrispettive annotazioni funzionali (.gff), estratte direttamente dalla banca dati NCBI.
A ciascun individuo nel [dataset](../dataset.tsv), le cui informazioni sono state reperite diettamente da NCBI, ci si riferisce ad essi come "specie NCBI".

## Download delle specie
Vengono scaricate le specie tramite uno specifico script
```bash
bash download_dataset.sh species_list_absent_in_GAGA.tsv
```


## Isoforma e cds
Dopo il download dei vari file, si è proceduto con il mantenimento dell'isoforma più lunga per ciascun gene identificato dall'annotazione e si è proceduto con l'estrazione delle cds.

Il tutto è stato eseguito mediante uno [script di snakemake](./longest.sh) in cui sono stati lanciati consecutivamente Vari programmi di agat.



```bash
snakemake -s longest.sh --cores 12 --use-conda
```


## Eliminazione degli pseudogeni
Tremite uno [script](./pseudogene_finde_eliminate.sh) si è proceduto con l'individuazione e l'eminiazione di tutte quelle sequenze identificabili come pseudogeni
```bash
bash ../01_NCBI_downloads/pseudogene_finde_eliminate.sh 
```


## Modifica deglie header
Anche per le specie NCBI è stato modificato ogni singolo header appartenente alle singole isoforme per ciascun gene.
```bash
for f in *.faa; do ID=$(basename -s .faa "$f"); sed -i -E "/^>/ s/>(.[^ ]+) gene=gene-(.[^ ]+) .+/>${ID}\|\2/" "$f"; done
```
