# Correzione e normalizzazione del p-value
Dopo che l'analisi della pgls è stata eseguita, si è proseguito con la normalizzazione del p-value per ognuno dei 296 file prodotti.
La normalizzazione (controllare se si può dire così) è stata eseguita con il successivo comando facente uso del metodo BH.
In questo modo nella colonna del p-value sono stati soprascritti i risultati dell'aggiustamento FDR ai precedenti.
```bash
for i in MA*; do 
    Rscript -e "
        df <- read.csv('$i/pgls5_castes.csv')
        df\$FDR <- p.adjust(df\$P_value, method='BH')
        
        # Reorder so FDR is column 4 (after the first 3 columns)
        rest <- setdiff(names(df)[-c(1:3)], 'FDR')
        df <- df[, c(names(df)[1:3], 'FDR', rest)]
        
        write.csv(df, '$i/pgls5_castes_adj.csv', row.names=FALSE)
    "
done
```
Per facilitare la visione e le successive operazioni tutti i file risultanti della pgls legate al fenotipo delle caste sono state unite in un unico grande file mediante l'uso del successivo comando:
```bash
awk -F',' '
BEGIN { OFS="\t"; header_written=0 }
{
    # Extract motif folder name from path (e.g. "MA0001.1/pgls5_castes_adj.csv" -> "MA0001.1")
    motif = FILENAME
    sub(/\/[^\/]+$/, "", motif)
    sub(/^.*\//, "", motif)

    # Remove all double quotes and convert commas to tabs
    gsub(/"/, "")
    gsub(/,/, OFS)

    # Print header only once from the first file, followed by data rows
    if (FNR == 1) {
        if (!header_written) {
            print $0, "motif"
            header_written = 1
        }
    } else {
        print $0, motif
    }
}' MA*/pgls5_castes_adj.csv > merged_pgls5_castes.tsv
```

# OGs significativi della PGLS
Al termine del lancio della pgls, il numero di motivi analizzati ammontava ad un valore estremamente grande con una altrettanto elevata varietà di risultati ottenuti. Per cercare di scremare l'elevato numero di OGs risultanti, dopo chiaramente aver normalizzato il parametro del p-value per ogni risultato di motivo, si è proceduto con il lancio di due differenti script mediante il seguente comando:

```bash
Rscript OG_sig_pvalue_only.R

Rscript OG_significativi_pgls5.R
```

Quello che è stato fatto dopo la nromalizzazione è stato selezionare tutti gli ortogruppi che presentassero un valore di p-value minore di 0.05 ed un valore di R^2 adj che veniva incrementato maggiormente di 0.15 ad ogni successiva esecuzione del programma (0.25, 0.40, 0.55, 0.70). Questo è necessario per restringere ulteriormente il campo dei possibili risultati, arrivando ad ottenere un numero lievemente ridotto di OG ma con una significatività estremamente elevata.

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

