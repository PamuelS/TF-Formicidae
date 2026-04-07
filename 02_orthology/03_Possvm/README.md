# Possvm 
Per poter eseguire il programma [Possvm](https://github.com/xgrau/possvm-orthology) in maniera corretta, è stato necessario creare e modificare alcuni file (prodotti da OrthoFinder e non). Quest oprogramma consente di eseguire una ulteriore scrematura dei geni ortologhi ottenuti in output dalla precedente analisi di OrthoFinder.

Per prima cosa si è proceduto con la modifica del file di "Resolved_gene_tree.txt" prodotto dall'analisi di OrthoFinder mediante una semplice sostituzione (dal momento che il programma prende in input solamente i file in formato .nwk)
```bash
sed -E 's/OG[0-9]+: //' Possvm_resolved_gene_tree.txt > Possvm_resolved_gene_tree.nwk
```

In seccondo luogo, per una migliore  riuscita dell'analisi di Possvm, si proceduto con la costruzione di un file che contenesse al programmad iinterpretare coerentemente gli abbreviativi che gli vengono forniti in input.
```bash
grep ">" *.fasta | sed 's/.*://; s/>//' | sed 's/^[^|]*|//' | awk '{id=$1; $1=""; sub(/^[ \t]+/, ""); print id "\t" $0}' > all_references.tsv
```

Una volta temrinate queste due operazioni, si posseggono tutti i file necessari per avviare l'analisi di Possvm, che è stata eseguita utilizzando il seguente codice:
```bash
python3 possvm.py -i Possvm_resolved_gene_tree.nwk  -spstree ../00_Orthofinder_analysis/OrthoFinder/Results_Mar30_1/Species_Tree/SpeciesTree_rooted.txt -r all_references.tsv -p possvm_ -split "|" -inflation 1.5 -o .
```
