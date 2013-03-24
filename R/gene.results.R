gene.results <-
function (object, cat.name = NULL, error = "none",
          print.it = TRUE, gene.names = NULL)
{
## gene.names: optional char vector of gene annotation
    if (!error %in% c("none", "FWER.Bonf", "FWER.Holm", "FDR.BH")) {
        cat(paste("WARNING: error = \"", error, "\" not available for local statistics\n",
            sep = ""))
        error <- "none"
    }
    C.names <- names(object@global.stat)
    if (is.null(cat.name))
        cat.name <- C.names[order(object@global.pval)][1]
    if (error != "none")
        error.p <- get(paste("error", error, sep = "."))(t(object@local.pval))
    keep <- as.matrix(object@C.mat[, C.names == cat.name])[,1]
    if (sum(is.na(object@local.pval)) > 0) {
        table <- cbind(Local.Stat = round(object@local.stat,
            3))[keep == 1, , drop = FALSE]
        table <- table[order(-table[, 1] * sign(table[, 1])),
            , drop = FALSE]
    }
    else if (error == "none") {
        table <- data.frame(Local.Stat = round(object@local.stat, 3),
            Emp.pvalue = sigfig(object@local.pval, 4), Temp = object@local.pval,
            stringsAsFactors = FALSE)[keep == 1, , drop = FALSE]
        table <- table[order(table[, 3], partial = -table[, 1] *
            sign(table[, 1])), , drop = FALSE][, 1:2, drop = FALSE]
    }
    else {
        table <- data.frame(Local.Stat = round(object@local.stat,3),
            Emp.pvalue = sigfig(object@local.pval, 4), Adj.pvalue = sigfig(error.p,4),
            Temp = object@local.pval, stringsAsFactors = FALSE)[keep == 1, , drop = FALSE]
        table <- table[order(table[, 4], partial = -table[, 1] *
            sign(table[, 1])), , drop = FALSE][, 1:3, drop = FALSE]
    }
    if (!is.null(gene.names)) {
      table <- data.frame(Gene.Names = gene.names[rownames(table)],
                          table, stringsAsFactors = FALSE)
    }
    if (print.it == TRUE) {
        cat("Category gene-specific results:\n")
        cat(paste("  Local:", object@local, "\n"))
        cat(paste("  Method:", object@method, "\n"))
        cat("\n")
        cat(paste(cat.name, "consists of", sum(keep), "gene features\n"))
        cat("\n")
        if (min(object@local.stat) < 0) {
            cat("Upregulated Genes\n")
            cat("-----------------\n")
        }
        if (max(table[, 'Local.Stat']) > 0)
            print(table[table[, 'Local.Stat'] > 0, , drop = FALSE])
        else cat("None in category\n")
        if (min(object@local.stat) < 0) {
            cat("\n")
            cat("Downregulated Genes\n")
            cat("-------------------\n")
            if (min(table[, 'Local.Stat']) <= 0)
                print(table[table[, 'Local.Stat'] <= 0, , drop = FALSE])
            else cat("None in category\n")
        }
    }
    else {
        table <- data.frame(Genes = dimnames(table)[[1]], table,
            stringsAsFactors = FALSE)
        rownames(table) <- NULL
        return(list(TableUp = as.list(table[table[, 'Local.Stat'] > 0, ,
            drop = FALSE]), TableDown = as.list(table[table[,
            'Local.Stat'] <= 0, , drop = FALSE])))
    }
}
