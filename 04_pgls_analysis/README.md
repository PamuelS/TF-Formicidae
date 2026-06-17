# PGLS analisi sui fenotipi
È stata eseguita la analisi pgls (Phylogenetic Generalized Least Square) per poter stimare concretamente la correlazione tra i motivi e i fenotipi raccolti dentro il dataset.
Nello specifico, si è cercato di correlare il punteggio dei motivi, sviluppato per ciascun ortogruppo associato ad ogni specie (ottenuto in [questo precedente passaggio](../03_tf_research/snakemake_motif_study.smk)), con lo tutti i disponibili stati variabili dei fenotipi investigati in questo studio (score ~ fenotipi).
Sono dunque state eseguite 8 differenti analisi pgls con lo scopo di correlare 6 distinti fenotipi:
1) 
