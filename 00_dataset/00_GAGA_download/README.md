# Specie GAGA
Vengono definite come "specie GAGA" tutti quegli individui presenti nel mio [dataset](../dataset.tsv) che sono stati inseriti anche nel dataset di GAGA.


## Download genomi 
sono stati scaricati i genomi dal [sito](https://sid.erda.dk/cgi-sid/ls.py?share_id=fU0yBp3NH5&current_dir=01_Genome_and_annotations&flags=f). Su questo sito sono state messe a disposizione tutti i dati relativi alle specie utilizzate nel dataset GAGA.

```bash
curl -L -o output ""https://sid.erda.dk/share_redirect/fU0yBp3NH5/01_Genome_and_annotations/01_Genome_and_annotations_allfiles.tar.gz
```
successivamente i file sono stati unzippati
```bash
tar -xzvf *tar.gz
```


## Modifica dei file 
una volta scaricati tutti i file relativi ad ogni specie nel dataset, si è proceduto con la modifica del nome della cartella associata alla specie e ad ognuno file di nostro interesse, che sono riassunti qui sotto:
- GAGA-0256_final_annotation_repfilt_addreannot_noparpse_representative.pep.fasta		(diventa: .faa)
- GAGA-0256_final_annotation_repfilt_addreannot_noparpse_representative_v3fixed.gff3	(diventa: _representative.gff3)
- GAGA-0256_SLR-superscaffolder_final_dupsrm_filt.repeats.gff				(diventa: _complete.gff)
- GAGA-0256_SLR-superscaffolder_final_dupsrm_filt.softMasked.fasta			(diventa: .fna)

è stato utilizzato il nome dei file gff e faa e fna per aggiungere l'abbreviativo della specie davanti al nome (per specie gaga):
```bash
for i in *; do cd "$i"; while read -r GAGA id; do mv "$GAGA"*_representative.pep.fasta "$id".faa; mv "$GAGA"*_representative_v3fixed.gff3 "$id"_representative.gff3; mv "$GAGA"*_filt.repeats.gff "$id"_complete.gff; mv "$GAGA"*_filt.softMasked.fasta "$id".fna; done < <(cut -f5,6 /DATASMALL/samuel.pederzini/TF-Formicidae/00_dataset/dataset.tsv | grep "GAGA" | tail -n+2); cd ..; done
```

è stato utilizzato il nome dei file gff e faa e fna per aggiungere l'abbreviativo della specie davanti al nome (per specie ncbi):
```bash
for i in *; do cd "$i"; while read -r ncbi id; do mv "$ncbi"*_representative.pep.fasta "$id".faa; mv "$ncbi"*_representative_v3fixed.gff3 "$id"_representative.gff3; mv "$ncbi"*_dupsrm_filt.repeats.gff "$id"_complete.gff; mv "$ncbi"*_dupsrm_filt.softMasked.fasta "$id".fna; done < <(cut -f5,6 /DATASMALL/samuel.pederzini/TF-Formicidae/00_dataset/dataset.tsv | grep "NCBI" | tail -n+2); cd ..; done
```

## Modifica degli header
Sono stati modificati anche gli header delle isoforme più lunghe di ogni gene, per ciascuna delle specie GAGA.
Si è quindi proceduto con la standardizzazione degli header ritrovati nei vari file .faa. 
Per riuscire ad associare l'abbreviativo GAGA al nostro abbreviativo personale, sono stati uniti due file mediante una colonna in comune
```bash
join <(cut -f1,4 list_species_GAGA.tsv) <(cut -f5,6 ../dataset.tsv | sort) > join_dataset.tsv
```

Questa standardizzazione consisteva nel posizionare l'abbreviativo utilizzato da GAGA all'inizio dell'header in modo tale che fosse più semplice identificarlo.


```bash
while read -r gaga abb; do [ -d "$abb" ] && { echo "Processing $abb..."; sed -i -E "/^>${gaga}/! s/^>(.*)_${gaga}_(.*)$/>${gaga}_\1_\2/" "$abb/$abb.faa"; } || echo "Errore: $abb non trovata"; done < <(cut -f2,3 /DATASMALL/samuel.pederzini/TF-Formicidae/00_dataset/00_GAGA_download/GAGA_vs_personal_ID.tsv | tr -d '\r')
```
dopo di che si proceduto con la sostituzione dell'abbreviativo GAGA con lìabbreviativo utilizzato da noi in questo studio 
```bash
while read -r gaga abb; do cd "$abb"; sed -i -E "/^>/ s/>${gaga}_/>${abb}\|/" "$abb".faa; cd ..; done < <(cut -f2,3 /DATASMALL/samuel.pederzini/TF-Formicidae/00_dataset/00_GAGA_download/GAGA_vs_personal_ID.tsv | tail -n+2)
```
