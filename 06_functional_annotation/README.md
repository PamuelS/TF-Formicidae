# Functional annotation
Per eseguire una annotazione funzionale degli ortogrippi si è inizialmente partiti con il raggruppamento di tutti gli ortogruppi che risultassero presenti nella analisi della pgls eseguita in precedenza.
Il suddetto raggruppamento di ortogruppi prenderà il nome di Gene-Universe al termine di tutto il procedimento di associazione dei Go term all'elenco di ortogruppi selezionato.

```bash
for i in ../05_pgls_analysis/01_pgls_results/MA*; do awk -F',' 'NR>1 && $9 ~ /success/ {gsub(/"/,"",$1); print $1}' "$i/pgls5_castes.csv"; done | sort -u > ortho_list.txt

```
Quindi sostanzialmente il gene-universe, che verrà ottenuto al termine di tutte le successive analisi, sarà costituito solo ed unicamente dagli OGs che hanno superato la precedente scrematura eseguita per l'analisi pgls.
Per questo procedimento è necessario individuare per ogni ortogruppo selezionato, la proteina più lunga ed associare ad essa la sequenza amminoacidica corrispondente. Questo procedimento viene svolto per ogni orogruppo preso in considerazione con il comando precedente e mediante lo script bash `longest_proteine.sh`.
```bash
bash longest_proteine_OG.sh
```
Il file risultante `longest_proteine_OGs.txt` conterrà una sola sequenza, associata alla specie di riferimento, per tutti gli ortogruppi elencati precedentemente.
Questo file viene poi sottoposto alle analisi di Diamond, per verificare a quale proteina quella sequenza possa essere associata, e ad InterProScan che identificherà tutti i possibili GOTerms riconosciuti in quella sequenza
Per avviare l'analisi eseguita da InterProScan è stato lanciato il seguente comando:
```bash
#per il background DISCO
interproscan.sh -i longest_protein_OGs.txt -goterms -pa -b missing -cpu 50 -applcommando

#per il background Ortho
interproscan.sh -i missing_longest_protein_ortho_OGs.txt -goterms -pa -b missing -cpu 50 -applcommando
```


Una volta ottenuto il file di annotazione funzionale generato direttamente da InterProScan, si può procedere all'eleiminazione di tutte quelle informazioni in eccesso fornite dal programma, per mostrare all'interno del file `go_back.tsv` solamente i GO che sono stati associati ad ogni singolo ortogruppo selezionato in precedenza.
```bash
awk -F'\t' '{
  gsub(/@.*/,"",$1); gsub(/\([^)]*\)/,"",$2); gsub(/\|/,",",$2);
  split($2, a, ",");
  for(i in a) if(a[i]!="") seen[$1,a[i]]=1
}
END {
  for(k in seen){
    split(k, b, SUBSEP)
    groups[b[1]] = (groups[b[1]] ? groups[b[1]] "," b[2] : b[2])
  }
  for(g in groups) print g "\t" groups[g]
}' <(cut -f1,14 samuel_motif_ants.tsv) | grep -v "-" > go_back.tsv

```

## Categorizzazione degli OGs significativi
Una volta ottenuto l'elenco di tutti gli ortogruppi che sono risultati significativi dai passaggi precedenti, si è proceduto con la categorizzazione di essi sulla base dl coefficente angolare riportato dall'analisi della pgls


```bash
Rscript extract_topgo_target.R
```
## Arricchimento funzionale
Al termine della selezione degli ortogruppi e della loro suddivisione dei due gruppi funzionali legati allo stato del fenotipo che meglio rappresenta il valore del coefficente angolare associato ad essi, si è proseguito con l'analisi dell'arrichimento funzionale dei GO terms:

```bash
for dir in 0*_enriched_*/0*_analysis; do
    echo "Running TopGO on: $dir"
    Rscript generalized_topGO.R "$dir"
done
```

Dal momento che i risultati potevano essere inseriti in differenti tipologie di casistich ebasate sulla tipologia di parametro scelto (p-value/FDR oppure bound specifico di R^2 adjusted) si è optato per la crezione di un unico file che consentisse di visualizzare tutte le informazioni in un singolo file eseguendo questo comando:
```bash
outfile="master_topGO_enrichment_results.tsv"

# 1. Grab header from the first available result file
first_file=$(find 0*_enriched_* -type f -name "topGOe_*" | head -n 1)

if [ -z "$first_file" ]; then
    echo "Error: No topGO output files found in subfolders!"
else
    # Write master header with metadata columns
    printf "R_squared\tAnalysis_Type\tTarget_Group\tSource_File\t" > "$outfile"
    head -n 1 "$first_file" >> "$outfile"

    # 2. Loop through all result files and merge content
    find 0*_enriched_* -type f -name "topGOe_*" | while read -r filepath; do
        filename=$(basename "$filepath")
        
        # Parse R^2 value from folder path (e.g., 0.25 from 01_enriched_0.25)
        r2=$(echo "$filepath" | sed -n 's/.*enriched_\([0-9.]*\).*/\1/p')
        
        # Parse analysis type
        case "$filepath" in
            *pvalue*) analysis="P-value" ;;
            *fdr*)    analysis="FDR" ;;
            *)        analysis="Unknown" ;;
        esac

        # Parse target group
        case "$filename" in
            *Polymorphic*) target="Polymorphic" ;;
            *Monomorphic*) target="Monomorphic" ;;
            *)             target="All" ;;
        esac

        # Prepend metadata columns and append rows (skipping individual file headers)
        tail -n +2 "$filepath" | awk -v r2="$r2" -v ana="$analysis" -v tar="$target" -v fn="$filename" \
            'BEGIN{FS="\t"; OFS="\t"} NF{print r2, ana, tar, fn, $0}' >> "$outfile"
    done

    echo "Done! Merged results saved to $outfile"
fi
```

## Associazione dell'OG al GO
Adesso si può eseguire l'associazione inversa per mappare quali sono gli ortogruppi associabili ai GO terms ottenuti dopo l'arricchimento. 
Il tutto viene eseguito partendo dal file di `go_back.tsv` che contiene tutti i possibili GO terms associabili a tutti gli ortogruppi selezionati per l'analisi.
```bash
Rscript extract_GO_genes.R
```


> Tutte queste istruzioni e comandi sono stati adoperati alla parte di dati risultanti dalla pipeline di DISCO. 
> Ricordarsi per ciò di applicarli nei minimi dettagli anche ai risultati della pipeline di ORTHO modificandone i path ed eventuali nomi


