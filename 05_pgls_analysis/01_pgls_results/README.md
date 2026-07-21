# OGs significativi della PGLS
Al termine del lancio della pgls, il numero di motivi analizzati ammontava ad un valore estremamente grande con una altrettanto elevata varietà di risultati ottenuti. Per cercare di scremare l'elevato numero di OGs risultanti, si è proceduto con il lancio di sue differenti script mediante il seguente comando:

```bash
Rscript OG_sig_pvalue_only.R

Rscript OG_significativi_pgls5.R
```

## OGs p-value significativi
Il primo script lanciato serve per selezionare solo ed unicamente tutti gli OGs che rispettassero il criterio di possedere un p-value che sia minore di 0.05. I risultati ottenuti dopo il lancio dello script mostano come il numero di ortogruppi totali a rispettare tale criterio sia passato a 240886 (da un valore iniziale di 2878598) avendo sempre una copertura totale di motivi pari a 296, overo tutti.

Per visualizzare meglio tale distribuzione, e per cercare di capire la direzione verso la quale un dererminato motivo puntava (se polimorfico oppure monomorfico) è stato eseguito anche il secondo script inerente a questo filtraggio

```bash
Rscript script_for_pgls_pvalue_dist.R
```

## OGs p-value && R^2 adj significativi
Successivamente è stato introdotto un livello ulteriore di selezione, rappresentato non più solo dal p-value minore di 0.05, ma anche dall'R^2 adjusted che doveva possedere un valore superiore a 0.25. In questa occasione, eseguendo sempre il secondo script come mostrato nel primo paragrafo, si ha avuto una drastica riduzione del numero di OGs selezionati arrivando a un totale di 1726 con una copertura in temrini di motivi paragonabile a 290 (6 motivi mancati dal totale).

Ad ogni modo come per il caso precedente, l'esecuzione di uno script apposito ha consentito di visualizzare accuratamente la distribuzione degli OGs all'interno dei motivi, permettendo anche di osservare la direzione vro la wuale un determinato motivo sta puntando.

```bash
Rscript script_pgls_total_sig_dist.R
```


> Tutto questo lavoro è stato svolto solo ed unicamente per la pgls legata all'analisi delle caste. Tutto il materiale riportato qui sopra potrà essere adoperato anche per le altre analisi pgls che dovranno essere eseguite per gli altri fenotipi.

La rappresentazione di tutti gli OG significativi per questi due parametri è stata riportata [qui](./OG_pvalue_Rsquaredadj_signif/pgls5_significant_OGs.tsv) sotto forma di tabella

