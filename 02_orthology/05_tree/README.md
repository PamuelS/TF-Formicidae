# Analisi filogenetica delle specie 
A partire da una topologia già nota (ovvero quella estratta dal papaer di GAGA menzionato in precedenza) vengono eseguite una seriue di operazioni in modo tale da aggiungere le specie NCBI non presenti nel dataset originario di GAGA ed eliminare tutte le specie categorizzate come Outgroup nel loro studio.

## Download degli alberi di GAGA
Per prima cosa sono stat iscaricati due alberi ottenuti dallo studio di GAGA
```bash
#file ottenuto da loro mediante allineamento dei codoni provenienti da 1000 geni
curl -L -o 1000genes ""https://sid.erda.dk/share_redirect/fU0yBp3NH5/05_Phylogeny/05b_Ortholog_trees/output_singlecopy_genesort1000_codon_nopartition_b1000_dnafixmodel_iqtree.treefile

#file ottenuto da loro mediante l'utilizzo di Astral
curl -L -o astral_database ""https://sid.erda.dk/share_redirect/fU0yBp3NH5/05_Phylogeny/05b_Ortholog_trees/singlecopy_codon_alltrees_bbcollapsed.astral.tree
```

## Modifica dei file
Una volta scaricati i file vengono tutti modificati in modo da sostituire il nome GAGA con l'abbreviativo scelto per questo lavoro. Inoltre vengono eliminati i bootstrap e i suporti ai nodi 
```bash
#eseguito su file 1000 geni
while read -r GAGA abb; do sed -i "s/${GAGA}/${abb}/g" 1000genes_modificato.nwk; done < <(cut -f5,6 /DATASMALL/samuel.pederzini/TF-Formicidae/00_dataset/dataset.tsv | tail -n+2)
sed -i -E 's/:[0-9.]+//g; s/[0-9./]+//g' 1000genes_modificato.nwk
```

A questo punto si è proceduto con l'inserimento delle specie NCBI assenti nel database GAGA. Il criterio di inserimento si è basato sull'eliminazione di tutte le politomie che si potevano creare durante la fase di inserimento delle specie. Le uniche politomie che sono state mantenute/create sono relative ai generi Lasius e Acromyrmex, dal momento che le esatte relazioni filogenetiche non erano ricreabili a causa di un cospiquo quantitativo di specie per genere.
L'inserimento delle specie ha tenuto conto, il più coerentemente possibile, delle informazioni filogenetiche reperite in bibliografia principalemnte legate ad analisi eseguita da altri studi sui rapporti filogenetici dei livelli tassonomi ci tribù.

Il file ottenuto al termine è `1000genes_modificato.nwk`

## Creazione dell'albero
Simultaneamente alla modifica della topologia preesistente, è stato lanciato il programma iqtree per costruire la filogenesi dei 175 campioni mediante l'utilizzo dei 200 ortogruppi `OG_tree.txt` che sono stati precedentemente allineati e trimmati.
```bash
iqtree -s ../04_trimmed/01_analysis_trimmed -m MFP+MERGE -B 1000 -T 4
```
