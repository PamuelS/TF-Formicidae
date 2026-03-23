# TF-Formicidae
Questa repository contiene tutti i file e i comandi che sono stati utilizzati per investigare la presenza di pattern regolatori differenziati nelle varie specie di formiche utilizzate in questa analisi.
Nello specifico quello che si è cercato corrisponde ad una eventuale espansione o contrazione dei motif (porzioni di genom,a di picole dimensioni alle qauli si attaccano i fattori trascrizionali) che possa essere collegata ad una variazione espressionale che porti alla codifica di un fenotipo differente, generalizzato per ogni specie di formica che possiede tale fenotipo.

I fenotipi investigari sono:
1) Quantità di regine all'interno del nido (Monoginia, Poliginia e Gamergate)
2) Caste (Monomorfica, Dimorfica e Polimorfica)
3) Parassitismo (Si/No)
   - Tipologia di parassitismo (Dulosi, Temporaneo, Inquilinismo e Lestobiosi)
4) Tipologia di ali della regina (Alata, Ergatoide e Brachiptera)
5) Nutrimento di emolinfa larvale (Si/No)
6) Bozzolo (Si/No)


## Workflow
1) Creazione del dataset e [download dei genomi e dei file .gff](./00_dataset)
2) Esecuzione del'[analisi BUSCO](./01_BUSCO)
3) Inferenza di [ortologia](./02_orthology)
