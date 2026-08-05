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

