# Functional annotation
Per eseguire una annotazione funzionale degli ortogrippi si è inizialmente partiti con il raggruppamento di tutti gli ortogruppi che risultassero presenti nella analisi della pgls eseguita in precedenza.
Il suddetto raggruppamento di ortogruppi prenderà il nome di Gene-Universe al termine di tutto il procedimento di associazione dei Go term all'elenco di ortogruppi selezionato.

```bash
for i in ../05_pgls_analysis/01_pgls_results/MA*; do awk -F',' 'NR>1 && $9 ~ /success/ {gsub(/"/,"",$1); print $1}' "$i/pgls5_castes.csv"; done | sort -u > ortho_list.txt

```
Quindi sostanzialmente il gene-universe è costituito solo ed unicamente dagli OGs che hanno superato la precedente scrematura eseguita per l'analisi pgls.
