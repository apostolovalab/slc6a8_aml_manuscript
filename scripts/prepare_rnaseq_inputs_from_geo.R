# Rebuild the per-run input files expected by */code/01-differential_expression.Rmd
# from the GEO count matrices (GEO: GSExxxxxx).
#
# Input  (per run): <run>/data/geo/RMB1_R<n>_raw_counts.tsv.gz   (from GEO)
#                   <run>/metadata/geo_samples.tsv                (in this repository)
# Output (per run): <run>/data/<ID>_results/<ID>-expression-matrix.tsv
#                   <run>/data/<ID>_results/<ID>-mapping-stats-reads.csv
#
# Run from the repository root:  Rscript scripts/prepare_rnaseq_inputs_from_geo.R

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# ==== Runs ====================================================================
runs <- tribble(
  ~project,                          ~geo_counts,                   ~run_id,
  "mrnaseq_rmb1_slc6a8_ko",          "RMB1_R1_raw_counts.tsv.gz",   "6GWQSF",
  "mrnaseq_rmb1_slc6a8_ko_run2",     "RMB1_R2_raw_counts.tsv.gz",   "KL7BFC",
  "mrnaseq_rmb1_slc6a8_ko_aza_ifng", "RMB1_R3_raw_counts.tsv.gz",   "7BD6RX"
)

# ==== Functions ===============================================================
read_samples <- function(project) {
  read_tsv(file.path(project, "metadata", "geo_samples.tsv"), show_col_types = FALSE) |>
    mutate(sample_num = as.integer(str_extract(sample_id, "(?<=_)[0-9]+$")),
           descriptor = str_remove(raw_file, "\\.fastq\\.gz$")) |>
    arrange(sample_num)
}

# Run 1 sample descriptors start at the sequencing index ("S1_Control_...")
format_descriptor <- function(descriptor, run_id) {
  if (run_id == "6GWQSF") str_remove(descriptor, paste0("^", run_id, "_[0-9]+_")) else descriptor
}

write_inputs <- function(project, geo_counts, run_id) {
  samples <- read_samples(project)
  counts  <- read_tsv(file.path(project, "data", "geo", geo_counts), show_col_types = FALSE)
  stopifnot(all(samples$library_name %in% colnames(counts)))

  expr <- counts |>
    select(gene_id, gene_name, gene_biotype, all_of(samples$library_name)) |>
    rename_with(~ paste0(samples$sample_id[match(.x, samples$library_name)], "_count"),
                all_of(samples$library_name))

  mapping <- tibble(SampleDescriptor = format_descriptor(samples$descriptor, run_id))

  out_dir <- file.path(project, "data", paste0(run_id, "_results"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(expr, file.path(out_dir, paste0(run_id, "-expression-matrix.tsv")))
  write_csv(mapping, file.path(out_dir, paste0(run_id, "-mapping-stats-reads.csv")))
  message(project, ": ", nrow(expr), " genes x ", nrow(samples), " samples")
}

# ==== Run =====================================================================
purrr::pwalk(runs, write_inputs)
