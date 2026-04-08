# Possvm 
Per poter eseguire il programma [Possvm](https://github.com/xgrau/possvm-orthology) in maniera corretta, è stato necessario creare e modificare alcuni file (prodotti da OrthoFinder e non). Quest oprogramma consente di eseguire una ulteriore scrematura dei geni ortologhi ottenuti in output dalla precedente analisi di OrthoFinder.

Per prima cosa si è proceduto con la modifica del file di `Resolved_gene_tree.txt` prodotto dall'analisi di OrthoFinder mediante una semplice sostituzione (dal momento che il programma prende in input solamente i file in formato .nwk)
```bash
sed -E 's/OG[0-9]+: //' Possvm_resolved_gene_tree.txt > Possvm_resolved_gene_tree.nwk
```

In seccondo luogo, per una migliore  riuscita dell'analisi di Possvm, si proceduto con la costruzione di un file che consentisse al programma di interpretare correttamente gli abbreviativi che gli vengono forniti in input e i nomi delle sequenze associate, il tutto contenuto in un unico file `all_references.tsv`.
```bash
grep ">" *.fasta | sed 's/.*://; s/>//' | sed 's/^[^|]*|//' | awk '{id=$1; $1=""; sub(/^[ \t]+/, ""); print id "\t" $0}' > all_references.tsv
```

In fine quello che è stato fatto è di spezzettare in molteplici subfile (circa 26000) l'albero `Possvm_resolved_gene_tree.nwk` dal momento che il programma non riusciva a leggere un file di tali dimensioni. Quindi si è proceduto con la separazione del file in molteplici sottoparti
```bash
split -l 1 -d -a 5 Possvm_resolved_gene_tree.nwk singoli_alberi/tree_

for f in tree_*; do mv "$f" "$f".nwk; done
```

Una volta temrinate queste due operazioni, si posseggono tutti i file necessari per avviare l'analisi di Possvm, che è stata eseguita utilizzando il seguente codice. Il problema è stato aggirato spezzettando in tante righe il file "Possvm_resolved_gene_tree.nwk" originario, eseguendo così l'analisi su ogni segmento del file.
```bash
ulimit -s unlimited; find singoli_alberi/ -name "*.nwk" | xargs -I {} -P 20 sh -c 'prefisso=$(basename {} .nwk); python3 possvm.py -i {} -spstree ../00_Orthofinder_analysis/OrthoFinder/Results_Mar30_1/Species_Tree/SpeciesTree_rooted.txt -o possvm_results/ -p "${prefisso}_" -split "|" -inflation 1.5'
```

I risulatati ottenuti per questa analisi, sono successivamente stati congiunti e accorpati, dal moemnto che il file originario `Possvm_resolved_gene_tree.nwk` era di dimensioni spropositate da riempiere immediatamente la memoria del computer
```bash
awk 'FNR==1 && NR!=1{next;}{print}' possvm_results/*.csv > RISULTATO_FINALE_ORTOLOGHI.csv
```
