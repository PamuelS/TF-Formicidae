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


# Annotazione funzionale dei risultati Orthofinder
La costruzione del background per l'annotazione funzionale dei risultati Ortho, è stata eseguita a partire dagli ortogruppi che mancavano nella creazione del background per i risultati DISCO. Di questa serie di OGs mancanti, è stata estratta la sequenza più lunga direttamente dagli ortogruppi che erano stati originati dall'analisi di OrthoFinder.

Quindi inizialmente è stata creata una lista di Ogs che possero presenti dentro l'analisi della pgls di Ortho, ma che mancassero dentro il file `go_back.tsv` relativo all'annotazione funzionale dei risultati DISCO.
```bash
while read -r og; do     [[ -z "$og" ]] && continue
    fa_file=$(find "$ORTHO_DIR" -name "${og}.fa*" | head -n 1)     if [[ -n "$fa_file" && -f "$fa_file" ]]; then         awk -v og="$og" '
            /^>/ {
                if (seq && length(seq) > maxlen) {
                    maxlen = length(seq)
                    maxhdr = hdr
                    maxseq = seq
                }
                hdr = $0
                sub(/^>/, "", hdr)
                # Take only the first word of the header if there are spaces
                split(hdr, a, " ")
                hdr = a[1]
                seq = ""
                next
            }
            { 
                # Remove gaps if any exist
                gsub(/-/, "", $0)
                seq = seq $0 
            }
            END {
                if (length(seq) > maxlen) {
                    maxhdr = hdr
                    maxseq = seq
                }
                if (maxhdr != "") {
                    printf ">%s@%s\n%s\n", og, maxhdr, maxseq
                }
            }
        ' "$fa_file" >> "$OUTPUT_FILE";     else         echo "Warning: Sequence file for $og not found in $ORTHO_DIR" >&2;     fi; done < missing_ortho_OGs_for_annotation.txt
```

> `missing_ortho_OGs_for_annotation.txt` è semplicemente l'elenco di ortogrppi mancanti ai quali bisogna associare la sequenza più luga ritrovata in tale ortogruppo


Dopo di che si è proceduto con l'associazione della proteina più lunga per ciascun ortogruppo selezionato 
```bash
ORTHO_DIR="/DATASMALL/samuel.pederzini/TF-Formicidae/02_orthology/00_Orthofinder_analysis/OrthoFinder/Results_Mar30_1/Orthogroup_Sequences"
OUTPUT_FILE="missing_longest_protein_ortho_OGs.txt"

> "$OUTPUT_FILE"

while read -r og; do
    [[ -z "$og" ]] && continue

    # Locate the matching FASTA file
    fa_file=$(find "$ORTHO_DIR" -name "${og}.fa*" | head -n 1)

    if [[ -n "$fa_file" && -f "$fa_file" ]]; then
        awk -v og="$og" '
            /^>/ {
                if (seq && length(seq) > maxlen) {
                    maxlen = length(seq)
                    maxhdr = hdr
                    maxseq = seq
                }
                hdr = $0
                sub(/^>/, "", hdr)
                split(hdr, a, " ")
                hdr = a[1]
                seq = ""
                next
            }
            { 
                gsub(/-/, "", $0)
                seq = seq $0 
            }
            END {
                if (length(seq) > maxlen) {
                    maxhdr = hdr
                    maxseq = seq
                }
                if (maxhdr != "") {
                    printf ">%s@%s\n%s\n", og, maxhdr, maxseq
                }
            }
        ' "$fa_file" >> "$OUTPUT_FILE"
    fi
done < missing_ortho_OGs_for_annotation.txt
```

Ed infine è avvenuta la ricostruzione completa del file gene universe per i risultati orthofinder, accorpando le informazioni provenienti dai file `go_back.tsv` (dove sono stati utilizzati gli ortogruppi presenti sia in DISCO  che in OrthoFinder) e dal file `go_back_missing.tsv` (che è stato creato dall'annotazione di InterProScan degli OGs che non erano presenti dentro il file gene universe di DISCO)

```bash
awk -F'\t' -v OFS='\t' '
$1 ~ /^OG/ {
    # Extract base OG name
    og = $1
    sub(/[:_].*/, "", og)
    
    # Split GO terms and keep uniques
    n = split($2, terms, ",")
    for (i=1; i<=n; i++) {
        term = terms[i]
        gsub(/^[ \t]+|[ \t]+$/, "", term)
        
        if (term != "") {
            if (!(og SUBSEP term in seen)) {
                seen[og SUBSEP term] = 1
                if (og in merged) {
                    merged[og] = merged[og] ", " term
                } else {
                    merged[og] = term
                }
            }
        }
    }
}
END {
    # Print the final combined list
    for (og in merged) {
        print og, merged[og]
    }
}' go_back.tsv go_back_missing.tsv | sort > go_back_ortho.tsv
```
