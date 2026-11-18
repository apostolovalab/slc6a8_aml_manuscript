h1 <- header_lvl(1) 

h2 <- header_lvl(2)

h3 <- header_lvl(3)

default_dt <- function(x, opts = NULL) {
  def_opts <- list(
    dom = "Blrtip",
    deferRender = FALSE,
    autoWidth = FALSE,
    buttons = list(
      I("colvis"),
      "csv",
      "excel"
    )
  )
  
  DT::datatable(
    x,
    rownames = FALSE,
    escape = FALSE,
    filter = "top",
    extensions = c("Buttons"),
    class = "compact stripe",
    width = "100%",
    options = c(def_opts, opts)
  )
}

contrast_to_group <- function(fit, contrast_matrix = NULL, cntr = NULL) {
  if (is.null(fit$contrasts)) {
    if (is.null(cntr) && is.null(contrast_matrix)) {
      stop("You specified an non-contrast fit object without specifying the
         contrast matrix and contrast of interest.")
    }
  }
  
  if (is.null(fit$contrasts)) {
    used_contrasts <- names(contrast_matrix[contrast_matrix[, cntr] != 0, cntr])
  } else {
    used_contrasts <- rownames(fit$contrasts)[fit$contrasts != 0]
  }
  
  contrast_names <- names(attributes(fit$design)$contrasts)
  
  group_list <- lapply(contrast_names, function(x) {
    if (any(grepl(x, used_contrasts))) {
      sub(x, "", used_contrasts)
    } else {
      used_contrasts
    }
  })
  
  names(group_list) <- contrast_names
  group_list[!vapply(group_list, is.null, FALSE)]
  return(group_list)
}

perform_contrast_fit_edgeR <- function(fit, cntr, contrast_matrix, filter_genes = NULL) {
  if (!is.null(filter_genes)) {
    fit <- fit[filter_genes(fit), ]
  }
  
  lrt <- glmQLFTest(fit, contrast = contrast_matrix[, cntr])
  lrt$contrasts <- contrast_matrix[, cntr, drop = FALSE]
  lrt$contrastCoef <- cntr
  lrt$fitPhenoData <- contrast_to_group(fit, contrast_matrix, cntr)
  return(lrt)
}

add_camera_stats <- function(fit, eset = NULL) {
  if (class(fit) == "MArrayLM" && is.null(eset)) {
    stop("Please supply expression data in the form of ExpressionSet using eset argument.")
  }
  
  if (class(fit) == "MArrayLM" && all(is.na(exprs(eset)))) {
    return(list())
  }
  
  if (class(fit) == "DGEGLM" && is.null(eset)) {
    stop("Please supply expression data in the form of DGEList using eset argument.")
  }
  
  fit$gseaStatistics <- cameraStats(
    y = eset,
    design = fit$design,
    contrast = fit$contrasts
  )
  return(fit)
}

perform_gsea <- function(fit, gene_set, min = 10, max = 200) {
  if (is.null(fit$table$logFC)) {
    warning("Data doesn't contain gene LogFC values, LogFC data will not be available for pathways.")
  }
  
  if (is.null(fit$gseaStatistics)) {
    return(data.frame())
  } else {
    indx <- gene_set %>%
      ids2indices(identifiers = rownames(fit$gseaStatistics$U))
    
    indx <- indx[map_int(indx, length) %>% between(min, max)]
    
    res <- cameraStats2Table(
      cameraStat = fit$gseaStatistics,
      index = indx,
      geneLog2FC = fit$table$logFC
    )
    return(res)
  }
}

expression_set_correlation <- function(eset, var_cut = 0.75, na = "drop") {
  dfCor <- as.data.frame(eset)
  
  if (na == "drop") {
    dfCor <- drop_na(dfCor)
  } else if (na == "replace_min") {
    dfCor[is.na(dfCor)] <- min(dfCor, na.rm = TRUE)
  } else {
    stop(
      c(
        "Incorrect value for how to deal with NA values. ",
        "Has to be one of 'drop', 'replace_min'."
      )
    )
  }
  
  iqrs <- apply(dfCor, 1, IQR)
  sel <- iqrs > quantile(iqrs, var_cut)
  
  return(cor(dfCor[sel, ]))
}

