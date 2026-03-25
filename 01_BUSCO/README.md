# Analisi di BUSCO
È stata eseguita una analisi BUSCO per tutte le specie NCBI, dal moemnto che non erano state rese disponibili online.
Non è stata eseguita l'analisi BUSCO sulle speci GAGA dal momento che i dati erano già stati resi disponibili da altri studi e ciascuna specie presentava volori dell'analisi accettabili.

BUSCO è stato lanciato tramite uno script snakemake per ottimizzare al meglio le risorse computazionali
```bash
snakemake -s ../snakemake_busco.smk --cores 12 --use-conda 
```

I risultati dell'analisi di BUSCO sono osservabili in questo [file](./00_busco_NCBI/all_busco_NCBI_species)
