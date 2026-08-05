# Analisi della pgls per OrhtoFinder
Le operazioni svolte sono le medesime, proprio per consetire un migliore confronto tra i risultati ottenuti dallo studio dei motivi con gli ortogruppi post DISCO e quelli ottenuti da OrthoFinder

## Applicazione correzione FDR

```bash
for i in MA*; do
    Rscript -e "
        df <- read.csv('$i/pgls5_Ortho_castes.csv')
        df\$FDR <- p.adjust(df\$P_value, method='BH')
        
        # Reorder so FDR is column 4 (after the first 3 columns)
        rest <- setdiff(names(df)[-c(1:3)], 'FDR')
        df <- df[, c(names(df)[1:3], 'FDR', rest)]
        
        write.csv(df, '$i/pgls5_Ortho_castes_adj.csv', row.names=FALSE)
    "
done
```

## Accorpamento dei risultati in un file

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
}' MA*/pgls5_Ortho_castes_adj.csv > merged_pgls5_Ortho_castes.tsv
```
