# Creazione del dataset
In questa cartella sono presenti tutti i genomi e lel loro relative annotazioni delle spèecie contenute all'interno del [dataset](./dataset.tsv). Sono state create due distinte categorie basate sulla presenza o assenza delle specie all'interno del dataset di GAGA (che fa riferimento a questo articolo: Adaptive radiation and social evolution of ants).

La determinazione dell'abbreviativo associato a ciascuna specie è stato scelto seguendo il criterio per il quale venivano selezionate le prime tre lettere del genere e le prime tre lettere della specie, per poi unirle assieme. Per le specie che sono rappresentate da più individui dentro il dataset o con una identificazione a livello di genere (sp.), si è optato per una modifica del criterio per l'abbreviazione della specie. Infatti con specie "doppie": un individuo eredità l'abbreviativo classico, mentre l'altro avrà l'abbreviativo caratterizzato da quattro lettere del genere e due lettere della specie.



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
