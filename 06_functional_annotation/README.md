# Functional annotation
Per eseguire una annotazione funzionale degli ortogrippi si è inizialmente partiti con il raggruppamento di tutti gli ortogruppi che risultassero presenti nella analisi della pgls eseguita in precedenza.
Il suddetto raggruppamento predne il nome di Gene-Universe ovvero l'insieme complessivo di ogni ortogruppo/gene utilizzato per questa analisi che fungerà da base per i successivi procedimenti.

```bash
for i in ../05_pgls_analysis/01_pgls_results/MA*; do awk -F',' 'NR>1 && $9 ~ /success/ {gsub(/"/,"",$1); print $1}' "$i/pgls5_castes.csv"; done | sort -u > gene_universe.txt

```
Quindi sostanzialmente il gene-universe è costituito solo ed unicamente dagli OGs che hanno superato la precedente scrematura eseguita per l'analisi pgls.
