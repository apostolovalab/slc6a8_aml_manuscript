# Reference files (not tracked)

Placeholders for public reference files used by the code. Download them and place them as below.

| Path | Used by | What to use |
|---|---|---|
| `references/MSigDb_2025.1/msigdb_v2025.1.Hs_gene_sets.rds` | BeatAML2, TCGA, figures | MSigDB v2025.1 human gene sets (https://www.gsea-msigdb.org), saved as an R data frame with one row per gene set and columns `collection_name`, `geneset_name`, `ensembl` (list column of Ensembl gene IDs). |
| `references/MSigDb_2025.1/msigdb_v2025.1.Mm_gene_sets.rds` | RMB-1 RNA-seq, figures | Same for the MSigDB v2025.1 mouse collections. |
| `atacseq_rmb1_slc6a8ko_aza_ifng/data/ref/GRCm39.primary_assembly.genome.fa.gz` | ATAC-seq (nf-core) | GENCODE mouse release M35, GRCm39 primary assembly genome (https://www.gencodegenes.org/mouse/). |
| `atacseq_rmb1_slc6a8ko_aza_ifng/data/ref/gencode.vM35.primary_assembly.annotation.gtf.gz` | ATAC-seq, ATAC figures | GENCODE M35 primary assembly annotation (GTF). |
| `atacseq_rmb1_slc6a8ko_aza_ifng/data/ref/mm39-blacklist.v2.liftover.bed` | ATAC-seq | ENCODE mm10 blacklist v2 (Amemiya et al. 2019) lifted over to mm39. |
| `atacseq_rmb1_slc6a8ko_aza_ifng/data/ref/jaspar2024_core_vertebrates_pwm.rds` | ATAC motif enrichment | Created by `03-motif_enrichment.Rmd` from the `JASPAR2024` Bioconductor package. |
| `beataml2_rnaseq/data/PanglaoDB_markers_27_Mar_2020.tsv` | BeatAML2 and TCGA figures | PanglaoDB marker table, 27 Mar 2020 release (https://panglaodb.se). |
| `beataml2_rnaseq/data/Cell_marker_Human.xlsx` | BeatAML2 and TCGA figures | CellMarker 2.0 human markers (http://bio-bigdata.hrbmu.edu.cn/CellMarker/). |
