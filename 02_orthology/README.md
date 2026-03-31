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
