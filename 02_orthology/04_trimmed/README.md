# Trimming delle sequenze
La stesa pulizia delle sequenze è avvenuta in un unico comando tramite l'utilizzo di snakemake, dato che sono stati eseguiti simultaneamente sia la funzionje di allineamento che di trimming mediante il lancio dello script `snakemake_alligned_trimmed.smk`.

Per dettagli guardare la [pagina precedente](../03_aligned)

Al termione del trimming delle sequenze, si è modificato la strutture degli header di ciascuno dei 200 ortogruppi selezionati per creare l'albero, in modo da creare una struttura di consenso da inciare ad iqtree

```bash
for i in *.faa; do sed -i -E "/^>/ s/\|.*$//" "$i"; done 
```
