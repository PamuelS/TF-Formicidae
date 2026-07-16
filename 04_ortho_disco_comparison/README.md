# Confronto tra ORTHOFINDER e DISCO
Come è possibile osservare nella [cartella precedente](../03_tf_research), è stata eseguita l'analisi dei motivi sia sugli ortogruppi che sono derivati dall'analisi di OrthoFinder che sui successivi ortogruppi generati dal lancio di DISCO.

Dal momento che questi due approcci di studio presentano criteri di generazione di ortogruppi estremamente differenziati, che sono riscontrabili proprio nel quantitativo di proteine/sequenze appartenenti ad ogni specie per il medesimo ortogruppo, questo ha fatto insorgere il dubbio che questo differente approccio potesse persino ripercuotersi nello studio dei motivi (come ad esempio una variazione nel punteggio di un motivo per un determinato ortogruppo). Perciò quello che è stato eseguito è un confronto diretto tra i risultati `totalscore_MAXXXX` ottenuti dall'analisi di OrthoFinder e quella DISCO sui medesimi ortogruppi per verificare una eventuale perdita di segnale derivata dall'eliminazione dei geni paraloghi all'interno degli ortogruppi generati da OrthoFinder.

Per eseguire con successo questo confronto sono state indette tre distinte categorie di ortogruppi basate sulla presenza o assenza di tale ortogruppo all'interno delle due distinte analisi.
- Ortho-Only ---> categori utilizzata per indicare gli OG che sono presenti unicamente dentro i file di ORTHOFINDER
- Split ---> categoria che indica gli OG presenti dentro i file di ORTHOFINDER e che hanno subito una ulteriore suddivisione (in varie forme alternative es. _00, _01, _02, ecc) dopo il lancio di DISCO
- Single-Complete ---> Categorai nella quale rientrano tutti gli OG che sono presenti nei file di entrambi i programmi e che non hannno subito ulteriori spacchettamenti da DISCO

## Cifre relative agli ortogruppi
L'esecuzione della analisi di Orthofinder ha riportato un quantitativo di ortofruppi paragonabile a 39226 (per vedere le statistiche in un formato maggiormente complessivo andare a guardare `Statistics_Overall.tsv`).
Gli ortogruppi che invece sono stati associati all'analisi dei motivi sono rappresentati al'incirca da 38000 OGs, con una minima riduzione del numero complessivo iniziale. Per rendere il confronto il più attendibile possibile, sono stati scartati da tutti i file di ORTHOFINDER gli ortogruppi che possedevano un quantitativo di specie inferiore a 90, passando così ad un ammonatare di 9789 OG. 

L'applicazione di DISCO ai risultati di ORTHOFINDER ha indubbiamente ridotto il numero di OG utilizzabili per l'analisi (dato che è stato utilizzato proprio il treshold di 90 specie per far si che DISCO prendesse un OG di ORTHOFINDER in considerazione), ma lo stesso numero è stato anche parzialmente ampliato anche grazie alla possibilità di DISCO di eseguire separazioni di OG. Infatti, la successiva esecuzione di DISCO sugli ortogruppi geneari da OrthoFinder ha prodotto ridotto il numero generale degli OG fino a 9783.
Di quei 9789 OG risultanti dopo il filtraggio per lo studio dei motvi, 9630 sono quelli che presentano una corrispondenza diretta con quelli DISCO e solamente 159 sono quelli che rimangono unicamente nei file ORTHOFINDER (I così detti Ortho-Only)

Sempre riferendosi agli OG di DISCO, quelli che ritroviamo sono proprio 9783 complessivi, entro i quali sono incliuse anche tutte le forme alternative generate dallo split che sono 84 e che sono riconducibili ad unicamente 72 OG (corrispondenti alla categoria Split)

## Anali del confronto
Mediante questo confronto sono state osservate due caratteristiche principali, ovvero il rapporto del punteggio complessivo di un ortogruppo tra DISCO e tra OERTHOFINDER ed inoltre il quantitativo di specie per un dato ortogruppo, che perdono il segnale con il motivo dopo l'esecuzione di DISCO.
Per il rapporto DISCO/ORTHO nel concreto viene eseguita la somma di tutti i punteggi ottenuti per il medesimo OG per ciascuno dei due file e viene successivamente eseguita la partizione di tale punteggio (qual'ora il risultato corrispondesse ad uno, significherebbe che non vi sono differenze per quel ortogruppo tra DISCO e ORTHOFINDER). Questo procedimento viene eseguito per ogni singolo OG ritrovato dentro al file di DISCO e viene riportato un singolo valore, corrispondente alla mediana di tutti i raporti, identificativo per un singolo motivo.
All'interno del grafico il boxplot corriponde al range del rapporto tra i due punteggi, mentre il violin plot posto in sottofondo mostra la frequenza di tale rapporto nei 296 motivi analizzati.
> Il rapporto DISCO/ORTHO è osservabile nel grafico [A]

L'altro aspetto analizzato per questo ocnfronto si è maggiormente concentrato sulla variazione della copertura delle specie che un motivo ha per un determinato OG. Nello specifico si è cercato di investigare se il passaggio da ORTHOFINDER a DISCO comportasse in alcuni ortogruppi la totale perdita di punteggio relativa ad una specie, questo perchè nel file ORTHOFINDER sono ammesse copie del medesimo gene associate ad una singola specie ma che vengono spacchettate e scomposte con il lancio di DISCO, causando una potenziale riduzione/perdita del punteggio complessivo di quel motivo per l'OG preso in considerazione.
Per calcolare questa differenza viene prima eseguita una stima in percentuale del quantitativo di specie con un punteggio maggiore di 0/NA per un OG (chiaramente sia per DISCO che per ORTHO) e viene ripetuto per ogni singolo OG disponbile. Successivamente viene calcolata una mediana complessiva dei valori ottenuti per tutti gli OG (sia per ORTHOFINDER che per DISCO), identificativa per un singolo motivo.
Infine viene calcolato il delta, corrispondente alla differenza della mediana tra DISCO e ORTHO per quel motivo.
> La copertura delle specie è osservabile nel grafico [B]
> 
> Un valore positivo di delta significa che l'esecuzione di DISCO non ha intaccato il quantitativo complessivo delle specie 
