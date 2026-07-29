# Correzione e aggiustamento del p-value
Dopo che l'analisi della pgls è stata eseguita, si è proseguito con un aggiustamento del p-value per ognuno dei 296 file prodotti.
L'aggiustamento è stato eseguito con il successivo comando facente uso del metodo BH.
In questo modo è stata aggiunta una colonna denominata FDR, affiancata a quella del p-value ottenuto dalla pgls analisi, utile per confrontare le differenze create dall'aggiustamento del parametro.
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

## Filtraggio degli OGs
Dopo questo procedimento si è proceduto con la selezione di tutti quegli OGs che rispecchiavano i parametri da noi imposti.
I criteri di selezione si basano su tre parametri fondamentali:
- soglia p-value
- soglia FDR
- soglia R^2 adjusted
I primi due parametri sono sempre rimast invariati, in modo che consentissero la selezione di OG che avessero sempre un parametro p-value o FDR minore di 0.05 (soglia di significaticità), mentre il terzo parametro è stato fatto variare, partendo con tutti i valori maggiori di 0.25 di R^2 adj, per poi passare a 0.40, 0.50 ed infine 0.70.

Quindi sostanzialmente sono state portate avanti due analisi parallele (p-value siginificativo da un lato e FDR significativo dall'altro) dove il parametro condiviso era il progressivo valore di R^2 adj che viene incrememtnato progressivamente. In questo modo si può valuatare lo stato delle due analisi confrontandole fra di loro
