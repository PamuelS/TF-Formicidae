Per lo studio dei motivi, è stato necessario estarre le informazioni da database online che contengono una serie di dati (come ad esempio le matrici di frequenza e di peso) necessari per approcciare uno studio di questo tipo. Tra i vari database esistenti, uno in particolare è stato anallizzato nel dettaglio, ovvero [TFLink](https://tflink.net/download/). Questo sito rappresenta una raccolta di qualsiasi database esistente che contenga informazioni relative ai Fattori Trascrizionali. Al suo interno sono contenuti tutti i 19 database esistenti, o quanto meno resi pubblici, ordinati secondo il cirterio di dimensione dove al primo posto si posiziona JASPAR 2024 CORE.


## Studio dei Motif
Per poter verificare la presenza dei motif all'interno dei 175 campioni utilizzati in questo studio, ci siamo affidati al database [JASPER 2026 CORE](https://jaspar.elixir.no/) aggiornato alla versione più recente 2026. Dopo una attenta analisi di tutti i databse presenti, si è optato per l'utilizzo dei dati provenienti univocamnete dal databse JASPAR 2026 CORE, dal momento che esso rappresenta la bancadati più ricca (per quanto riguarda i motif) tra tutti quelli resi disponibili su TFLink, ed inoltre si è evitato una potenziale ridondanza di motivi con nomenclature diverse provenienti da database differenti.
Nel sito JASPAR sono stati selezionati tutti e 357 i motif disponibili che sono stati associati alla specie *Drosophila melanogaster*. In seguito sono stati scaricati songolarmente tutti i motif mediante un singolo link ed è stato unzippato il file che li conteneva tutti
```bash
wget https://jaspar.elixir.no/temp/20260416114301_JASPAR2026_individual_matrices_2399676_jaspar.zip

unzip 20260416114301_JASPAR2026_individual_matrices_2399676_jaspar.zip 
```
