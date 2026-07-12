# Confronto tra ORTHOFINDER e DISCO
Come è possibile osservare nella [cartella precedente](../03_tf_research), è stata eseguita l'analisi dei motivi sia sugli ortogruppi che sono derivati dall'analisi di OrthoFinder che sui successivi ortogruppi generati dal lancio di DISCO.

Dal momento che questi due approcci di studio presentano criteri di generazione di ortogruppi estremamente differenziati, che sono riscontrabili proprio nel quantitativo di proteine/sequenze appartenenti ad ogni specie per il medesimo ortogruppo, questo ha fatto insorgere il dubbio che questo differente approccio potesse persino ripercuotersi nello studio dei motivi (come ad esempio una variazione nel punteggio di un motivo per un determinato ortogruppo). Perciò quello che è stato eseguito è un confronto diretto tra i risultati `totalscore_MAXXXX` ottenuti dall'analisi di OrthoFinder e quella DISCO sui medesimi ortogruppi per verificare una eventuale perdita di segnale derivata dall'eliminazione dei geni paraloghi all'interno degli ortogruppi generati da OrthoFinder.

Per eseguire con successo questo confronto sono state indette tre distinte categorie di ortogruppi basate sulla presenza o assenza di tale ortogruppo all'interno delle due distinte analisi.
- Ortho-Only ---> categori utilizzata per indicare gli OG che sono presenti unicamente dentro i file di ORTHOFINDER
- Split ---> categoria che indica gli OG presenti dentro i file di ORTHOFINDER e che hanno subito una ulteriore suddivisione (in varie forme alternative es. _00, _01, _02, ecc) dopo il lancio di DISCO
- Single-Complete ---> Categorai nella quale rientrano tutti gli OG che sono presenti nei file di entrambi i programmi e che non hannno subito ulteriori spacchettamenti da DISCO

