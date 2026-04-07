# DISCO
Prima di lanciare il comando di disco, sono stati modificati tutti i nomi delle specie del file Resolved_Gene_Trees.txt, che OrthoFinder ha ritoccato dopo l'esecuzione dell'analisi. Questo procedimento è necessario per la buona riuscita del programma DISCO.
```bash
sed -E 's/[A-Z][a-z]{5}_//g; s/\)n[0-9]*+/\)/g' Resolved_Gene_Trees.txt
```
A questo pèunto si può procedere con il lancio del programma DISCO, che eliminerà i geni paraloghi dagli ortogruppi formatisi dopo il lancio di OrthoFinder
```bash
while IFS='|' read -r OG tree; do python3 ../../../../01_DISCO/disco.py -i <(echo "$tree") -o ../../../01_DISCO/${OG/:/}.nwk -d "|" -m 90 --remove_in_paralogs --keep-labels --verbose >> ../../../01_DISCO/disco.log; done < Resolved_Gene_Trees.txt
```

In seguito all'operazione svolta dal programma DISCO, si è proceduto con il controllo dati eliminando per prima cosa i file risultanti come vuoti e spostandoli in una cartella specifica.
```bash
find . -size 0 -print > empty_disco.txt

find . -size 0 -delete
```

Sono stati salvati anche gli ortogruppi che al termine dell'operazione di DISCO presentano tutte le specie presenti nel dataset.
```bash
for i in *.nwk; do [[ $(grep -oE '\b[A-Z][a-z]{5}\b' "$i" | wc -l) -eq 175 ]] && echo "$i"; done > every_species
```
