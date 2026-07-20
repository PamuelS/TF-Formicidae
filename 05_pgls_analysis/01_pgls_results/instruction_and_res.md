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
