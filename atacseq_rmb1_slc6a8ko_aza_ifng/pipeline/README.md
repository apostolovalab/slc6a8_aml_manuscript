# ATAC-seq preprocessing

nf-core/atacseq 2.1.2 (Nextflow 24.10.4). Download the FASTQ files from GEO into `data/fastq/`, then from `atacseq_rmb1_slc6a8ko_aza_ifng/`:

```bash
nextflow run nf-core/atacseq -r 2.1.2 -profile singularity -params-file pipeline/params.yaml
```

`samplesheet_nfcore.csv` lists the two sequencing lanes of each library (24 libraries, 48 rows). If your Nextflow setup requires it, replace the FASTQ paths with absolute paths.
