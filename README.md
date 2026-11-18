# SLC6A8 in AML – analysis code

Code to reproduce the analyses and figures of:

> Apostolova P, Gremelspacher T, Sonrel A, *et al.* The creatine transporter SLC6A8 is an immune pressure-exposed metabolic vulnerability in acute myeloid leukemia. *Journal* (Year). DOI: TBD

Code DOI: [10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX) · Processed data DOI: [10.5281/zenodo.YYYYYYY](https://doi.org/10.5281/zenodo.YYYYYYY)

## Data

| Data | Source | Identifier |
|---|---|---|
| RMB-1 bulk RNA-seq (runs 1–3) and ATAC-seq: raw reads, count matrices, peaks, bigWigs | GEO | GSExxxxxx |
| RMB-1 proteomics | ProteomeXchange / MassIVE | PXD087485 / MSV000104035 |
| Processed data (DiffBind counts, CIBERSORTx fractions, ProteoFlux output, …) | Zenodo | 10.5281/zenodo.YYYYYYY |
| LC-MS metabolomics input tables | Massive | TBD |
| BeatAML2 (RNA-seq, clinical) | BeatAML2 data portal | Bottomly et al., *Cancer Cell* 2022 |
| TCGA-LAML (RNA-seq, clinical, mutations) | NCI GDC (downloaded by `tcga_slc6a8_rnaseq/code/00_tcga_download.Rmd`) | TCGA-LAML |
| AML single-cell atlas | Whittle et al., *eLife* 2026 | see `manuscript_figures/data/README.md` |
| Tabula Sapiens blood | CELLxGENE Discover | dataset `1bf55482-ad49-446a-87c6-f3442bebcef1` |

Each analysis folder has a `data/README.md` listing the files it expects and where to get them. Reference files (genome, annotation, gene sets) go in `references/` (see `references/README.md`).

## Repository structure

Each analysis is a self-contained folder; R Markdown files are run from their `code/` folder (relative paths `../data`, `../results`).

```
mrnaseq_rmb1_slc6a8_ko/            RMB-1 RNA-seq run 1 (SLC6A8 KO vs NTC; vehicle, T cells, AZA, AZA + T cells)
mrnaseq_rmb1_slc6a8_ko_run2/       RMB-1 RNA-seq run 2 (vehicle, AZA + T cells)
mrnaseq_rmb1_slc6a8_ko_aza_ifng/   RMB-1 RNA-seq run 3 (AZA + IFN-gamma)
atacseq_rmb1_slc6a8ko_aza_ifng/    RMB-1 ATAC-seq (vehicle, AZA, AZA + IFN-gamma)
beataml2_rnaseq/                   BeatAML2 cohort
tcga_slc6a8_rnaseq/                TCGA-LAML cohort
manuscript_figures/                manuscript figures (incl. metabolomics, proteomics, single-cell atlas)
references/                        reference files (not tracked)
scripts/                           helper to build RNA-seq inputs from GEO
```

## Running the analyses

Run the files in this order; each step writes to its own `results/` folder, which the next steps read.

1. **RMB-1 RNA-seq** (each run folder): download the GEO count matrices into `<run>/data/geo/`, then from the repository root run `Rscript scripts/prepare_rnaseq_inputs_from_geo.R`. Then render `01-differential_expression.Rmd`, then `02-tf_inference.Rmd`. 
2. **ATAC-seq**: place `dba_counted.rds` (Zenodo) in `atacseq_rmb1_slc6a8ko_aza_ifng/results/rds/diffbind/`, then render `01-differential-accessibility.Rmd` and `02-motif_enrichment.Rmd`. Locus plots in `manuscript_figures/code/atacseq_figures.Rmd` also need the filtered BAM files and genome FASTA from nf-core/atacseq 2.1.2, run on the GEO FASTQ files with `atacseq_rmb1_slc6a8ko_aza_ifng/pipeline/`.
3. **BeatAML2**: `1a` → `1b` → `2b` (CIBERSORTx input files) → `2a` → `1ab` → `3a`. CIBERSORTx fractions (Zenodo) go in `beataml2_rnaseq/results/tables_slc6a8_continuous/`.
4. **TCGA-LAML**: `00` → `1b` → `1a` (CIBERSORTx input files) → `2a`. CIBERSORTx fractions (Zenodo) go in `tcga_slc6a8_rnaseq/results/tables_slc6a8_continuous/`.
5. **Figures**: render the files in `manuscript_figures/code/`.

Example: `cd beataml2_rnaseq/code && Rscript -e 'rmarkdown::render("1a_slc6a8_continuous_DE.Rmd")'`

### CIBERSORTx

Cell-type fractions were estimated with the CIBERSORTx web portal (https://cibersortx.stanford.edu) from the mixture and single-cell reference files written by `beataml2_rnaseq/code/2b_celltype_deconvolution.Rmd` and `tcga_slc6a8_rnaseq/code/1a_celltype_deconvolution.Rmd` (Module 1: Create Signature Matrix; Module 2: Impute Cell Fractions). The resulting fractions are provided on Zenodo.

## Software

R 4.5.2 with Bioconductor 3.22. ATAC-seq preprocessing: nf-core/atacseq 2.1.2 (Nextflow 24.10.4).

## License

Code: MIT (see `LICENSE`).

## Contact

Anthony Sonrel – anthony.sonrel@unibas.ch
