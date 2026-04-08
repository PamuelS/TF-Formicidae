# Split DISCO
Partendo dai risultati ottenuti dal lancio di DISCO, si procede con l'esecuzione dello [script](./alternative_split_disco_output.sh).

```bash
bash ../02_DISCO_OG/alternative_split_disco_output.sh /DATASMALL/samuel.pederzini/TF-Formicidae/02_orthology/00_Orthofinder_analysis/OrthoFinder/Results_Mar30_1/Orthogroup_Sequences
```

Nello specifico questo script consente di prednere l'output prodotto direttamente da [DISCO](../01_DISCO) (ovvero gli alberi in formato .nwk) ed estrae i nomi delle sequenze che ritrova all'interno dei file, per poi associarci direttamente le sequenze amminoacidiche che sono state ottenute dall'analisi di [OrthoFinder](../00_Orthofinder_analysis)

Dal momento che i dati di partenza possedevano annotazioni differenziate (essendo che erano stati utilizzati differenti programmi di annotazione da differenti persone), ho dovuto miodificare il file `split_disco_output.sh` modificando la parte relativa al comando `grep` e nella parte terminale del `awk`. Tutte le modifiche sono state apportate nel file `alternative_split_disco_output.sh` 
