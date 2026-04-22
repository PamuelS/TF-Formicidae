Per lo studio dei motivi, è stato necessario estarre le informazioni da database online che contengono una serie di dati (come ad esempio le matrici di frequenza e di peso) necessari per approcciare uno studio di questo tipo. Tra i vari database esistenti, uno in particolare è stato anallizzato nel dettaglio, ovvero [TFLink](https://tflink.net/download/). Questo sito rappresenta una raccolta di qualsiasi database esistente che contenga informazioni relative ai Fattori Trascrizionali. Al suo interno sono contenuti tutti i 19 database esistenti, o quanto meno resi pubblici, ordinati secondo il cirterio di dimensione dove al primo posto si posiziona JASPAR 2024 CORE.


# Studio dei Motif
Per poter verificare la presenza dei motif all'interno dei 175 campioni utilizzati in questo studio, ci siamo affidati al database [JASPER 2026 CORE](https://jaspar.elixir.no/) aggiornato alla versione più recente 2026. Dopo una attenta analisi di tutti i databse presenti, si è optato per l'utilizzo dei dati provenienti univocamnete dal databse JASPAR 2026 CORE, dal momento che esso rappresenta la bancadati più ricca (per quanto riguarda i motif) tra tutti quelli resi disponibili su TFLink, ed inoltre si è evitato una potenziale ridondanza di motivi con nomenclature diverse provenienti da database differenti.
Nel sito JASPAR sono stati selezionati tutti e 357 i motif disponibili che sono stati associati alla specie *Drosophila melanogaster*. In seguito sono stati scaricati songolarmente tutti i motif mediante un singolo link ed è stato unzippato il file che li conteneva tutti
```bash
wget https://jaspar.elixir.no/temp/20260416114301_JASPAR2026_individual_matrices_2399676_jaspar.zip

unzip 20260416114301_JASPAR2026_individual_matrices_2399676_jaspar.zip 
```

Sono state scaricate tutte le versioni dei motivi disponibili per *D. melanogaster*, perciò persino vecchie versioni. Si è proceduto con lo spostamento delle versioni più arretrate dentro ad una ulteriore cartella.
```bash
# spostamento versione 1
for i in *.2.jaspar; do motif=$(basename "$i" .2.jaspar); if [ -f "$motif".1.jaspar ]; then echo "sposto la versione vecchia del moivo $i"; mv "$motif".1.jaspar older_version/; fi; done

# spostamento versione 2
for i in *.3.jaspar; do motif=$(basename "$i" .3.jaspar); if [ -f "$motif".2.jaspar ]; then echo "sposto la versione vecchia del moivo $i"; mv "$motif".2.jaspar older_version/; fi; done
```
Al termine di questa operazione i motif risultanti corrispondono a 297

## Utilizzo di PWM Scan
I file ottenuti nel precedente passaggio, sono stati scaricati nel formato .JASPAR, corrispondente quindi ad una PFM (Position Frequency Matrices), che corrisponde proprio alla frequenza con la quale una determinata base viene ritrovata all'interno di un motivo in differenti campioni. Per poter utilizzare questo dato, è necessario convertitlo in una ulteriore forma di matrice probabilistica (che può essere PWM oppure PSSM).

Per eseguire la conversione da un modello di matrice ad un altro, è stato adoperato il programa [PWM Scan](https://sourceforge.net/projects/pwmscan/reviews/), che innanzitutto è stato scaricato in locale e poi trasferito sul server:
```bash
scp /home/STUDENTI/samuel.pederzini/Downloads/pwmscan.1.1.9.tar.gz STUDENTI^samuel.pederzini@137.204.142.237:/DATASMALL/samuel.pederzini/TF-Formicidae/03_tf_research
```
