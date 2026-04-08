# Split DISCO
Partendo dai risultati ottenuti dal lancio di DISCO, si procede con l'esecuzione dello [script](./split_disco_output).

```bash
bash ../02_DISCO_OG/split_disco_output.sh /DATASMALL/samuel.pederzini/TF-Formicidae/02_orthology/00_OrthoFinder_analysis/OrthoFinder/Results_Mar30_1/Orthogroup_sequences
```

Nello specifico questo script consente di prednere l'output prodotto direttamente da [DISCO](../01_DISCO) (ovvero gli alberi in formato .nwk) ed estrae i nomi delle sequenze che ritrova all'interno dei file, per poi associarci direttamente le sequenze amminoacidiche che sono state ottenute dall'analisi di [OrthoFinder](../00_Orthofinder_analysis)
