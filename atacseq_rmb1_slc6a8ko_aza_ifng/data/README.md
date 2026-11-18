# Input data – RMB-1 ATAC-seq

| File | Source |
|---|---|
| `../results/rds/diffbind/dba_counted.rds` | Zenodo (DiffBind object with the consensus peaks and raw fragment counts per sample). With this file, `02-differential-accessibility.Rmd` starts from the counts and does not need the BAM files. |
| `fastq/*.fastq.gz` | GEO GSExxxxxx (only needed to rerun nf-core/atacseq, see `../pipeline/`) |
| `ref/` | Genome, annotation and blacklist, see `references/README.md` |

nf-core/atacseq outputs (`../results/bowtie2/merged_library/*.mLb.clN.sorted.bam`, `../results/genome/`) are needed only to count reads from scratch (`params$recount = TRUE`) and for the locus plots in `manuscript_figures/code/atacseq_figures.Rmd`.
