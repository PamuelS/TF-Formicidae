# Orthology Inference
In questa cartella sono riportate tuttee le informazioni e i dati realtivi all'analisi di inferenza di ortologia per tutte le 175 specie di formiche.

Per prima cosa è stato lanciato il programma Orthofinder su tutti i proteomi delle formiche (sia GAGA che NCBI)
```bash
orthofinder -f ./whole_proteome -t 25
```
a causa dell'interruzione improvvisa dell'analisi, terminata dopo l'esecuzione del blast fra tutte le specie, si è proseguito con l'utilizzo di una ulteriore flag nel comando di orthofinder che ha consentito di riprendere l'analisi da dove si è interrotta, prendendo in input i file blast creati con l'analisi precedente.
```bash
orthofinder -b whole_proteome/OrthoFinder/Reults_Mar25/ -t 25
```

## Ulteriore modifica degli header
A causa di una mancata modifica degli header associati ad alune delle specie presenti nel dataset, è stata eseguita una modifica ulteriore e successiva degli header in tutti i file prodotti dall'analisi di orthofinder, mediante il seguente comando
```bash
while read -r gaga abb; do grep -rl "$gaga" . --exclude-dir=WorkingDirectory | xargs sed -i -E "s/${gaga}_?/${abb}\|/g"; done < <(cut -f 2,3 ../../../../00_dataset/00_GAGA_download/GAGA_vs_personal_ID.tsv | tail -n+2)
```
è stato eseguito un comando per verificare che non fossero state apportate modifiche persino alle sequenze amminoacidiche.
```bash
for i  in *; do grep -l "|" <(grep -v ">" "$i"); done
```
ed inoltre si è verificato che il quantitativo di ">" corrsipondesse al quantitativo di "|" mediante un ulteriore comando
```bash
for i in *; do pipe=$(grep -c "|" "$i"); great=$(grep -c ">" "$i"); if  [ "$pipe" -ne "$great" ]; then echo "$i"; fi; done
```
Una volta individuati qui file che possedevano ancora un header non standardizzato secondo i criteri scelti per questo studio, si è proceduto con la modifica a mano degli header
```bash
# per vedere quale header non possedeva l'annotazione corretta
diff -y <(grep ">" OG0000579.fa | sort | sed 's/>//') <(grep "|" OG0000579.fa | sort | sed 's/>//') | egrep "<|>"

# per associare l'header corretto alla specie di appartenenza
grep "XM_070293761" *.faa
```

