makeFitObject <- function(fits) { 
    
    fitObjList <- list()
    for (fn in names(fits)) {
        fitObj <- list()
        fit <- fits[[fn]]

        if (is(fit, "MArrayLM")) {
            require(limma)

            fitObj$samples <- rownames(fit$design)

            coef <- colnames(fit$coefficients)
            i <- grep("Intercept", coef)
            if (length(i)) {
                coef <- coef[-i]
            }
            if (length(coef) != 1)
                stop('Expect a single contrast coefficient per fit')
            fitObj$contrastCoef <- coef

            if (!is.null(fit$fitPhenoData))
                fitObj$fitPhenoData <- fit$fitPhenoData

            fitObj$meanExpression <- fit$Amean

            tt <- topTable(fit, coef=fitObj$contrastCoef, number=nrow(fit))
            fitObj$topOrder <- rownames(tt)
            tt <- tt[names(fitObj$meanExpression),]
            se <- sqrt(fit$s2.post) * fit$stdev.unscaled[,fitObj$contrastCoef]
            fitObj$statistics <- data.frame(logFC=tt[,'logFC',drop=FALSE],
                                            stdErr=se,
                                            P.Value=tt[,'P.Value',drop=FALSE],
                                            adj.P.Val=tt[,'adj.P.Val',drop=FALSE])

            fitObj$gseaStatistics <- zscoreT(fit$t[,coef],fit$df.total)
            fitObj$design <- fit$design
            fitObj$contrasts <- fit$contrasts

            
        } else if (is(fit, "DGELRT")) {
            require(edgeR)

            fitObj$samples <- rownames(fit$design)

            if (!is.null(fit$contrasts)) {
                coef <- colnames(fit$contrasts)
                i <- grep("Intercept", coef)
                if (length(coef) != 1)
                    stop('Expect a single contrast coefficient per fit')
                fitObj$contrastCoef <- coef
            }
            
            if (!is.null(fit$fitPhenoData))
                fitObj$fitPhenoData <- fit$fitPhenoData
            
            fitObj$meanExpression <- fit$AveLogCPM
            names(fitObj$meanExpression) <- rownames(fit$genes)

            tt <- topTags(fit, n=nrow(fit))$table
            fitObj$topOrder <- rownames(tt)
            tt <- tt[names(fitObj$meanExpression),]
            fitObj$statistics <- data.frame(logFC=tt[,'logFC', drop=FALSE],
                                            P.Value=tt[,'PValue'],
                                            adj.P.Val=tt[,'FDR'])
            fitObj$gseaStatistics <- fit$gseaStatistics$Stat
            fitObj$design <- fit$design
            fitObj$contrasts <- fit$contrasts

        } else {
            cl <- class(fit)
            stop(sprintf("\"%s\" type fit object is currently not supported!", cl))
        }
        fitObjList[[fn]] <- fitObj
    }
    fitObjList
}

