# PGLS analisi sui fenotipi
È stata eseguita la analisi pgls (Phylogenetic Generalized Least Square) per poter stimare concretamente la correlazione tra i motivi e i fenotipi raccolti dentro il dataset.
Nello specifico, si è cercato di correlare il punteggio dei motivi, sviluppato per ciascun ortogruppo associato ad ogni specie (ottenuto in [questo precedente passaggio](../03_tf_research/snakemake_motif_study.smk)), con lo tutti i disponibili stati variabili dei fenotipi investigati in questo studio (score ~ fenotipi).

Sono dunque state eseguite 8 differenti analisi pgls con lo scopo di correlare 6 distinti fenotipi:
1) pgls1 -> parassitismo (assente/presente)
2) pgls2 -> emolinfa larvale (assente/presente)
3) pgls3 -> bozzolo (assente/presente)
4) pgls4 -> ali della regina
  - pgls4a (alata/ergatoide)
  - pgls4b (alata/brachiptera)
6) pgls5 -> caste (monomorfiche/polimorfiche)
7) pgls6 -> numero di regine
  - pgls6a (monoginia/poliginia)
  - pgls6b (monoginia/facoltative/poliginia)

L'intera struttura dello script R lasi può trovare [qui](./pgls_phenotypes.R), mentre il codice snakemake utilizzato per gestire un elevato quantitativo di dati è reperibile in [questo](./snakefile_pgls_analysis.smk) punto.
