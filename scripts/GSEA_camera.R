## GSEA using camera 
require(limma)
require(edgeR)
## 
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param y 
##' @param design 
##' @param contrast 
##' @param weights 
##' @param allow.neg.cor 
##' @param trend.var 
##' @param fit 
##' @param ... 
##' @return 
##' @author Robert Ivanek
cameraStats <- function(y, design=NULL, contrast=ncol(design),
                        weights=NULL, trend.var=FALSE, use.ranks = FALSE,
                        inter.gene.cor = 0.01, fit=NULL, ...) {
  ## if y is "DGEList" then run edge R code
  if (is(y, "DGEList")) {
    y <- edgeR:::.zscoreDGE(y = y, design = design, contrast = contrast)
  }
  
  ## taken from the limma code
  y <- getEAWP(y)
  G <- nrow(y$exprs)
  n <- ncol(y$exprs)
  ID <- rownames(y)
  if (G < 3) 
    stop("Two few genes in dataset: need at least 3")
  if (is.null(design)) 
    design <- y$design
  if (is.null(design)) 
    stop("design matrix not specified")
  else {
    design <- as.matrix(design)
    if (mode(design) != "numeric") 
      stop("design must be a numeric matrix")
  }
  if (nrow(design) != n) 
    stop("row dimension of design matrix must match column dimension of data")
  p <- ncol(design)
  df.residual <- n - p
  if (df.residual < 1) 
    stop("No residual df: cannot compute t-tests")
  if (is.null(weights)) 
    weights <- y$weights
  fixed.cor <- !(is.na(inter.gene.cor) || is.null(inter.gene.cor))
  if (fixed.cor) {
    if (use.ranks) 
      df.camera <- Inf
    else df.camera <- G - 2
  } else {
    df.camera <- min(df.residual, G - 2)
  }
  y <- y$exprs
  if (!is.null(weights)) {
    if (any(weights <= 0)) 
      stop("weights must be positive")
    if (length(weights) == n) {
      sw <- sqrt(weights)
      y <- t(t(y) * sw)
      design <- design * sw
      weights <- NULL
    }
  }
  if (!is.null(weights)) {
    if (length(weights) == G) 
      weights <- matrix(weights, G, n)
    weights <- as.matrix(weights)
    if (any(dim(weights) != dim(y))) 
      stop("weights not conformal with y")
  }
  if (is.character(contrast)) {
    contrast <- which(contrast == colnames(design))
    if (length(contrast) == 0) 
      stop("coef ", contrast, " not found")
  }
  if (length(contrast) == 1) {
    j <- c((1:p)[-contrast], contrast)
    if (contrast < p) 
      design <- design[, j]
  }
  else {
    QR <- qr(contrast)
    design <- t(qr.qty(QR, t(design)))
    if (sign(QR$qr[1, 1] < 0)) 
      design[, 1] <- -design[, 1]
    design <- design[, c(2:p, 1)]
  }
  if (is.null(weights)) {
    QR <- qr(design)
    if (QR$rank < p) 
      stop("design matrix is not of full rank")
    effects <- qr.qty(QR, t(y))
    unscaledt <- effects[p, ]
    if (QR$qr[p, p] < 0) 
      unscaledt <- -unscaledt
  }
  else {
    effects <- matrix(0, n, G)
    colnames(effects) <- ID
    unscaledt <- rep.int(0, G)
    names(unscaledt) <- ID
    sw <- sqrt(weights)
    yw <- y * sw
    for (g in 1:G) {
      xw <- design * sw[g, ]
      QR <- qr(xw)
      if (QR$rank < p) 
        stop("weighted design matrix not of full rank for gene ", 
             g)
      effects[, g] <- qr.qty(QR, yw[g, ])
      unscaledt[g] <- effects[p, g]
      if (QR$qr[p, p] < 0) 
        unscaledt[g] <- -unscaledt[g]
    }
  }
  U <- effects[-(1:p), , drop = FALSE]
  sigma2 <- colMeans(U^2)
  U <- t(U)/sqrt(pmax(sigma2, 1e-08))
  if (trend.var) 
    A <- rowMeans(y)
  else A <- NULL
  sv <- squeezeVar(sigma2, df = df.residual, covariate = A)
  modt <- unscaledt/sqrt(sv$var.post)
  if (!is.null(fit) || is(fit, "MArrayLM")) {
    ## moderated T statistics including mixed effect
    Stat <- zscoreT(fit$t[,contrast], df = fit$df.total)
  } else {
    if (use.ranks) 
      Stat <- modt
    else {
      df.total <- min(df.residual + sv$df.prior, G * df.residual)
      Stat <- zscoreT(modt, df = df.total, approx = TRUE)
    }
  }
  ## 
  return(list(Stat=Stat, U=U, df.camera=df.camera))
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param cameraStat 
##' @param index 
##' @param use.ranks 
##' @param allow.neg.cor 
##' @param sort 
##' @param ... 
##' @return 
##' @author Robert Ivanek
cameraStats2Table <- function(cameraStat, index, use.ranks=FALSE,
                              allow.neg.cor=FALSE, inter.gene.cor = 0.01,
                              sort=TRUE, adjust.method="BH", geneLog2FC=NULL, ...) {
  if (!is.list(index)) 
    index <- list(set1 = index)
  nsets <- length(index)
  if (nsets == 0) 
    stop("index is empty")
  ## extract the components from the list
  Stat <- cameraStat$Stat
  U <- cameraStat$U
  df.camera <- cameraStat$df.camera
  G <- nrow(U)
  n <- ncol(U)
  ID <- rownames(U)
  ##
  fixed.cor <- !(is.na(inter.gene.cor) || is.null(inter.gene.cor))
  ##
  meanStat <- mean(Stat)
  varStat <- var(Stat)
  tab <- matrix(0, nsets, 5)
  rownames(tab) <- names(index)
  colnames(tab) <- c("NGenes", "Correlation", "Down", "Up", "TwoSided")
  for (i in 1:nsets) {
    iset <- index[[i]]
    if (is.character(iset)) 
      iset <- which(ID %in% iset)
    StatInSet <- Stat[iset]
    m <- length(StatInSet)
    m2 <- G - m
    if (fixed.cor) {
      correlation <- inter.gene.cor
      vif <- 1 + (m - 1) * correlation
    }
    else {
      if (m > 1) {
        Uset <- U[iset, , drop = FALSE]
        vif <- m * mean(colMeans(Uset)^2)
        correlation <- (vif - 1)/(m - 1)
      }
      else {
        vif <- 1
        correlation <- NA
      }
    }
    tab[i, 1] <- m
    tab[i, 2] <- correlation
    if (use.ranks) {
      if (!allow.neg.cor) 
        correlation <- max(0, correlation)
      tab[i, 3:4] <- rankSumTestWithCorrelation(iset, statistics = Stat, 
                                                correlation = correlation, df = df.camera)
    }
    else {
      if (!allow.neg.cor) 
        vif <- max(1, vif)
      meanStatInSet <- mean(StatInSet)
      delta <- G/m2 * (meanStatInSet - meanStat)
      varStatPooled <- ((G - 1) * varStat - delta^2 * m * 
                          m2/G)/(G - 2)
      two.sample.t <- delta/sqrt(varStatPooled * (vif/m + 
                                                    1/m2))
      tab[i, 3] <- pt(two.sample.t, df = df.camera)
      tab[i, 4] <- pt(two.sample.t, df = df.camera, lower.tail = FALSE)
    }
  }
  tab[, 5] <- 2 * pmin(tab[, 3], tab[, 4])
  tab <- data.frame(tab, stringsAsFactors = FALSE)
  Direction <- rep.int("Up", nsets)
  Direction[tab$Down < tab$Up] <- "Down"
  tab$Direction <- Direction
  ## add average abslotute log2 FC for every geneset
  if (!is.null(geneLog2FC)) {
    absLog2FC <- sapply(index, function(x) mean(abs(geneLog2FC[x]), na.rm=T))
    tab$absLog2FC <- absLog2FC[rownames(tab)]
  }
  tab$PValue <- tab$TwoSided
  tab$Down <- tab$Up <- tab$TwoSided <- NULL
  ## if (fixed.cor) 
  ##     tab$Correlation <- NULL
  if (nsets > 1) 
    tab$FDR <- p.adjust(tab$PValue, method = adjust.method)
  if (sort && nsets > 1) {
    o <- order(tab$PValue)
    tab <- tab[o, ]
  }
  tab
}
