# Functional annotation
Per eseguire una annotazione funzionale degli ortogrippi si è inizialmente partiti con il raggruppamento di tutti gli ortogruppi che risultassero presenti nella analisi della pgls eseguita in precedenza.
Il suddetto raggruppamento di ortogruppi prenderà il nome di Gene-Universe al termine di tutto il procedimento di associazione dei Go term all'elenco di ortogruppi selezionato.

```bash
for i in ../05_pgls_analysis/01_pgls_results/MA*; do awk -F',' 'NR>1 && $9 ~ /success/ {gsub(/"/,"",$1); print $1}' "$i/pgls5_castes.csv"; done | sort -u > ortho_list.txt

```
Quindi sostanzialmente il gene-universe, che verrà ottenuto al termine di tutte le successive analisi, sarà costituito solo ed unicamente dagli OGs che hanno superato la precedente scrematura eseguita per l'analisi pgls.
Per questo procedimento è necessario individuare per ogni ortogruppo selezionato, la proteina più lunga ed associare ad essa la sequenza amminoacidica corrispondente. Questo procedimento viene svolto per ogni orogruppo preso in considerazione con il comando precedente e mediante lo script bash `longest_proteine.sh`.
```bash
bash longest_proteine.sh
```
Il file risultante `longest_proteine.tsv` conterrà una sola sequenza, associata alla specie di riferimento, per tutti gli ortogruppi elencati precedentemente.
Questo file viene poi sottoposto alle analisi di Diamond, per verificare a quale proteina quella sequenza possa essere associata, e ad InterProScan che identificherà tutti i possibili GOTerms riconosciuti in quella sequenza

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

