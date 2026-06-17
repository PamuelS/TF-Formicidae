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

### Accorgimenti
- L'albero che è stato utilizzato è l'albero ottenuto con l'analisi di Maximum Likelihood eseguito fornendo la topologia dell'albero ottenuta dal lavoro di GAGA e con l'aggiunta delle 12 "specie NCBI" (ottenuto mediante la scelta di 500 geni DISCO).
- Il dataset è stato parzialmente modificato in modo che potesse restituire un quantitativo di stati alternativi dei fenotipi riconducibile ad un massimo 3/4 e qual'ora fosse stato possibile lo stato veniva ridotto ad una forma di presenza o assenza.
- Essendo che l'analisi pgls si basava sui motivi associati agli ortogruppi post DISCO, si è presupposto che il punteggio risultante fosse razzionalmente distribuito entro un range di valori non troppo elevati e con una bassa presenza di outlier. Perciò i valori sono stati successivamente normalizzati secondo il prcedimento statistico min-max entro un intervallo che varia da 0 a 1 e perciò non sono stati gestiti eventuali osservazioni discostanti dalla norma mediante la statistica MAD.
- La pgls è stata eseguita su ciascun ortogruppo che fosse quantomeno associato una singola volta con almeno una specie, il tutto eseguit per ciascun motivo.
- Sono stati imposti quattro differenti "bound", ovvero dei limiti prefissati di Lambda di Pagel (λ) che è stata calcolata successivamente mediante Maximum Likelihood. I bound imposti sono: (1e-05, 1), (1e-03, 1), (1e-01, 1), (1, 1).
- Sono stati inseriti dei limiti di accettazione, legati al quantitativo di specie che presentano quel fenotipo, per far si che l'analisi di uno specifico ortogruppo potesse essere preso in considerazione.
