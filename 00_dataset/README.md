# Creazione del dataset
In questa cartella sono presenti tutti i genomi e lel loro relative annotazioni delle spèecie contenute all'interno del dataset. Sono state create due distinte categorie basate sulla presenza o assenza delle specie all'interno del dataset di GAGA (che fa riferimento a questo articolo: Adaptive radiation and social evolution of ants).

Le specie appartenenti alla categoria "GAGA" sono già presenti nel dataset GAGA e sono state semplicemente scaricate da [qui](https://sid.erda.dk/cgi-sid/ls.py?share_id=fU0yBp3NH5&current_dir=01_Genome_and_annotations&flags=f). Per quanto riguarda le specie non presenti in tale papaer, esse sono state scaricate direttamente da NCBI, tramite l'uso di uno script apposito (download_dataset.sh).

```bash
bash download_dataset.sh .
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
