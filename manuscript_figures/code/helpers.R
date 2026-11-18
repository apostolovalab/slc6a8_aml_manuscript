summarise_mat <- function(m, label) {
  observed <- m[!is.na(m)]
  data.frame(
    Stage          = label,
    Features       = nrow(m),
    Samples        = ncol(m),
    Total_Values   = length(m),
    Observed       = sum(!is.na(m)),
    Missing        = sum(is.na(m)),
    Pct_Missing    = round(100 * mean(is.na(m)), 2),
    Mean           = round(mean(observed), 4),
    Median         = round(median(observed), 4),
    SD             = round(sd(observed), 4),
    Min            = round(min(observed), 4),
    Max            = round(max(observed), 4),
    Features_w_NA  = sum(apply(m, 1, anyNA)),
    Samples_w_NA   = sum(apply(m, 2, anyNA))
  )
}

make_heatmap_metab <- function(mat_name, dat_list, hits_core, hits_any) {  
  
  mat_raw  <- dat_list[[mat_name]]
  grp_raw  <- as.character(attr(mat_raw, "group"))
  grp      <- sub("^[^_]+_", "", grp_raw)
  
  all_hits <- unique(c(hits_core, hits_any[[mat_name]]))
  rows_use <- all_hits[all_hits %in% rownames(mat_raw)]
  
  mat_sub  <- mat_raw[rows_use, , drop = FALSE]
  mat_z    <- t(scale(t(as.matrix(mat_sub))))
  
  row_split    <- factor(ifelse(rows_use %in% hits_core, "Core", "Any"),
                         levels = c("Core", "Any"))
  row_fontface <- ifelse(rows_use %in% hits_core, "bold", "plain")
  
  grp_levels <- unique(grp)
  grp_cols   <- cs$SampleGroup
  
  col_anno <- HeatmapAnnotation(
    Group = anno_block(
      gp        = gpar(fill = grp_cols[grp_levels], col = NA),
      labels    = grp_levels,
      labels_gp = gpar(fontsize = 8, col = "black")
    ),
    simple_anno_size     = unit(5, "mm"),
    show_annotation_name = FALSE,
    show_legend          = FALSE
  )
  
  col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#D6604D"))
  
  Heatmap(
    mat_z,
    name                  = "Z-score",
    col                   = col_fun,
    
    column_title          = NULL,
    
    top_annotation        = col_anno,
    
    column_split          = factor(grp, levels = grp_levels),
    cluster_column_slices = FALSE,
    column_gap            = unit(1, "mm"),
    
    row_split             = row_split,
    row_title             = NULL,
    row_gap               = unit(3, "mm"),
    
    show_row_names        = TRUE,
    row_names_side        = "left",
    row_names_gp          = gpar(fontsize = 9, fontface = row_fontface),
    show_row_dend         = FALSE,
    cluster_rows          = TRUE,
    cluster_row_slices    = FALSE,
    
    cluster_columns       = FALSE,
    show_column_names     = FALSE,
    
    width  = unit(4 * ncol(mat_z), "mm"),
    height = unit(4 * nrow(mat_z) + 20, "mm")
  )
}

