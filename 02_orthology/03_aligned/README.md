# Allineamento delle sequenze
In questa cartella vengono svolte tutte le operazioni necessarie ad eseguire gli allineamenti delle sequenze dei geni ortologhi single-copy ottenuti dai passaggi precedenti.
Per raggiungere un numero di ortogruppi pari a 200, sono stati utilizzati tutti i `single-copy` ottenuti dall'[analisi di OrthoFinder](../00_Orthofinder_analysis/Statistics_Overall.tsv) e per raggiungere il numero di ortogruppi prestabilito sono stati scelti gli ortogruppi originatisi dopo il lancio delcomando `alternative_split_disco_outgroup.sh`

I mancanti 71 ortogruppi `single-copy` sono stati scelti randomicamente dai risultatti di DISCO
```bash
less every_species_OGs.txt | shuf -n 71 > species_tree.txt
```

> il file `every_species_OGs.txt` indica quanti sono gli ortogruppi che risultanti come single-copy al termine di `alternative_split_disco_outgroup.sh`

Per verificare se gli ortogruppi presi randomicamente non si ripetessero con i `single-copy` di OrthoFinder, è stato eseguito il successivo comando:
```bash
sed -E 's/_00//' species_tree_OG.txt | sort | uniq -c | less
```
