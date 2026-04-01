# DISCO
Prima di lanciare il comando di disco, sono stati modificati tutti i nomi delle specie del file Resolved_Gene_Trees.txt, che OrthoFinder ha ritoccato dopo l'esecuzione dell'analisi. Questo procedimento è necessario per la buona riuscita del programma DISCO.
```bash
sed -E 's/[A-Z][a-z]{5}_//g; s/\)n[0-9]*+/\)/g' Resolved_Gene_Trees.txt
```
