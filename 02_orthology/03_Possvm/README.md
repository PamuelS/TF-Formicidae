# Possvm 
Per poter eseguire il programma [Possvm](https://github.com/xgrau/possvm-orthology) in maniera corretta, è stato necessario creare e modificare alcuni file (prodotti da OrthoFinder e non).

Per prima cosa si è proceduto con la modifica del file di "Resolved_gene_tree.txt" prodotto dall'analisi di OrthoFinder mediante una semplice sostituzione (dal momento che il programma prende in input solamente i file in formato .nwk)
```bash
sed -E 's/OG[0-9]+: //' Possvm_resolved_gene_tree.txt > Possvm_resolved_gene_tree.nwk
```
