# Input data – RMB-1 RNA-seq run3

| File | Source |
|---|---|
| `geo/RMB1_R3_raw_counts.tsv.gz` | GEO GSExxxxxx (supplementary file of the RNA-seq series) |
| `7BD6RX_results/7BD6RX-expression-matrix.tsv`, `7BD6RX_results/7BD6RX-mapping-stats-reads.csv` | Created from the GEO file by `scripts/prepare_rnaseq_inputs_from_geo.R` (sample table: `../metadata/geo_samples.tsv`) |
| `collectri_mouse.rds` | CollecTRI mouse regulons from OmnipathR (code in the `eval=FALSE` chunk of the TF inference report; needs internet access) |
| `Stat1_chipatlas.tsv` | ChIP-Atlas "Target Genes" table for mouse Stat1 (https://chip-atlas.org/target_genes), columns `Target_genes` and `Stat1|Average` |
| `stat1_regulons_signor.csv` | STAT1 targets from SIGNOR (code in the `eval=FALSE` chunk of the TF inference report; needs internet access) |

Gene sets: `references/MSigDb_2025.1/` (see `references/README.md`).
