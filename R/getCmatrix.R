getCmatrix <-
function (keyword.list = NULL, gene.list = NULL,
          present.genes = NULL, GO.ont = NULL,
          min.size = 2, max.size = Inf, by.gene = FALSE,
          gene.names =  NULL, prune = FALSE,
          as.matrix = FALSE, ...)
{
##    ...: Allows arguments from versions 1.0 and 2.0 to be ignored
## keyword.list, gene.list: Objects of class 'list' containing category assigment by
##        column or by row respectively
## present.genes: Optional character vector of row labels to match C.mat to X.mat
## GO.ont: Optional character string, "CC","BP", or "MF" to specify ontology
## min.size, max.size: Optional numeric thresholds for gene-set size
## by.gene: Logical as to whether multiple probesets to a single gene should be downweighted
## gene.names: Character vector  of gene names required for by.gene = TRUE
## prune: Logical to remove redundant categories
## as.matrix: Logical to return an object of class 'matrix' instead of 'SparseM'

require(GO.db); require(SparseM)
  if(is.null(present.genes))   if(is.null(keyword.list)){
      if(prod(present.genes %in% names(gene.list)) == 0)
        cat("WARNING: Some present.genes do not appear in gene.list\n")
      if(sum(present.genes %in% names(gene.list)) == 0)
         stop("No present.genes matched the names of gene.list",call. = FALSE)
      present.genes <- names(gene.list)
  } else present.genes <- sort(unique(unlist(keyword.list)))
  if(by.gene){
     rowvec <- table(gene.names,exclude=NULL)
     rowvec <- 1/rowvec[match(gene.names,names(rowvec))]
     rowvec[is.na(names(rowvec))] <- 1
  } else rowvec <- rep(1,length(present.genes))

if(is.null(keyword.list)){
    C.names <- sort(unique(unlist(gene.list)))
    genes.match = match(present.genes, names(gene.list))
    present.list <- gene.list[genes.match]

    ra = 1; ja = 1; ia = 1; raja.index = 1
    cat("  Genes completed:")
    num <- length(present.list)
     for(g in 1:num){
       if (g %in% floor(seq(num/5,num,num/5))) cat(paste(ceiling(g/num*100),"% ",sep=""))
       gene.cats = present.list[[g]]
       if(!is.null(gene.cats)){
       if(all(is.na(gene.cats))) gene.cats = NA else {
         gene.cats <- unique(gene.cats[!is.na(gene.cats)]) }
         num.cats = length(gene.cats)
       if(!all(is.na(gene.cats))){
          cat.matches = match(gene.cats, C.names)
          last.index = raja.index + num.cats - 1
          ra[raja.index:last.index] = rep(rowvec[g], num.cats)
          ja[raja.index:last.index] = cat.matches
          ia[g] = raja.index
          raja.index = raja.index + num.cats
       } else ia[g] = raja.index } else ia[g] = raja.index
    }
    ia[g+1] = raja.index
    dimension = c(num, length(C.names))
    C.matrix.csr = new("matrix.csr", ra = ra, ja = as.integer(ja),
       ia = as.integer(ia), dimension = as.integer(dimension))
    rowvec = as.matrix.csr(matrix(rowvec, nrow = 1))
    size = as.matrix(rowvec %*% C.matrix.csr)
    C.matrix.csr <- C.matrix.csr[, (size >= min.size) & (size <= max.size)]
    C.names = C.names[(size >= min.size) & (size <= max.size)]
  } else {
    if(!is.null(GO.ont)){
      keep <- names(keyword.list) %in% names(as.list(GOTERM))
      keep2 <- sapply(mget(names(keyword.list)[keep],GOTERM),Ontology) %in% GO.ont
      keyword.list <- keyword.list[keep][keep2]
    }
    if(is.null(present.genes)) present.genes <- sort(unique(unlist(keyword.list)))
    keep <- rep(0,length(keyword.list))
    for (i in 1:length(keyword.list)) {
      keep[i] <- ((rowvec %*% (present.genes %in% keyword.list[[i]])) >= min.size) &
                 ((rowvec %*% (present.genes %in% keyword.list[[i]])) <= max.size)
    }
    keyword.list <- keyword.list[keep==1]
    C.names <- names(keyword.list)

    ra = 1;ja = 1; ia = 1; raja.index = 1
    cat("  Categories completed:")
    num <- length(keyword.list)
    for(g in 1:num){
       if (g %in% floor(seq(num/5,num,num/5))) cat(paste(ceiling(g/num*100),"% ",sep=""))
       gene.cats = unique(keyword.list[[g]][keyword.list[[g]] %in% present.genes])
       num.cats = length(gene.cats)
       if(!is.null(gene.cats)){ if(!is.na(gene.cats[1])){
          cat.matches = match(gene.cats, present.genes)
          last.index = raja.index + num.cats - 1
          ra[raja.index:last.index] = rowvec[cat.matches]
          ja[raja.index:last.index] = cat.matches
          ia[g] = raja.index
          raja.index = raja.index + num.cats
       } else ia[g] = raja.index } else ia[g] = raja.index
    }
    ia[g+1] = raja.index
    dimension = c(length(C.names), length(present.genes))
    C.matrix.csr = new("matrix.csr", ra = ra, ja = as.integer(ja),
       ia = as.integer(ia), dimension = as.integer(dimension))
    C.matrix.csr = t(C.matrix.csr)
  }
  cat(paste("\n ",length(C.names),"catgories formed\n"))

  if(prune){
     drop <- apply(as.matrix(C.matrix.csr)>0,2,which)
     drop <- duplicated(sapply(drop,paste,collapse=","))
     C.matrix.csr = C.matrix.csr[,!drop]
     C.names <- C.names[!drop]
     cat(paste("  Pruned to",sum(!drop),"catgories \n"))
  }

  if(as.matrix) {
    C.mat <- as.matrix(C.matrix.csr)
    dimnames(C.mat) <- list(present.genes,C.names)
    return(C.mat)
  } else return(list(C.mat.csr=C.matrix.csr, row.names = present.genes, col.names = C.names))
}
