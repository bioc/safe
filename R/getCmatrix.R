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

#### 0) If rows not specified, build from lists; else check for errors
  if(is.null(present.genes)){
    if(is.null(keyword.list)) present.genes <- names(gene.list) else
      present.genes <- sort(unique(unlist(keyword.list)))
  } else {
    if(is.null(keyword.list)) {
      if(prod(present.genes %in% names(gene.list)) == 0)
         cat("WARNING: Some present.genes do not appear in gene.list\n")
      if(sum(present.genes %in% names(gene.list)) == 0)
         stop("No present.genes match the names of gene.list",call. = FALSE)
    } else if(sum(present.genes %in% unlist(keyword.list)) == 0)
            stop("No present.genes match the elements of keyword.list",call. = FALSE)
  }
#### 1) Create weights 
  if(by.gene){
     rowvec <- table(gene.names,exclude=NULL)
     rowvec <- 1/rowvec[match(gene.names,names(rowvec))]
     rowvec[is.na(names(rowvec))] <- 1
  } else rowvec <- rep(1,length(present.genes))
#### 2) Make C from gene.list
if(is.null(keyword.list)){
    C.names <- sort(unique(unlist(gene.list)))
    genes.match = match(present.genes, names(gene.list))
    present.list <- gene.list[genes.match]

    drop <- sapply(present.list,is.null)
    present.list[!drop] <- lapply(present.list[!drop],function(x) x[!is.na(x)])
    present.list <- lapply(present.list, unique)
    num.cats.vec <- sapply(present.list, length)
    ra <- unlist(mapply(rep, rowvec, num.cats.vec))
    ja <- unname(unlist(lapply(present.list, match, C.names)))
    ia <- unname(cumsum(c(1,num.cats.vec)))

    dimension = c(length(present.list), length(C.names))
    C.matrix.csr = new("matrix.csr", ra = ra, ja = as.integer(ja),
       ia = as.integer(ia), dimension = as.integer(dimension))
    rowvec = as.matrix.csr(matrix(rowvec, nrow = 1))
    size = as.matrix(rowvec %*% C.matrix.csr)
    C.matrix.csr <- C.matrix.csr[, (size >= min.size) & (size <= max.size)]
    C.names = C.names[(size >= min.size) & (size <= max.size)]
} else {
#### 3) Make C from keyword.list
    if(!is.null(GO.ont)){
      keep <- names(keyword.list) %in% names(as.list(GOTERM))
      keep2 <- sapply(mget(names(keyword.list)[keep],GOTERM),Ontology) %in% GO.ont
      keyword.list <- keyword.list[keep][keep2]
    }
    
    keep <- rep(0,length(keyword.list))
    for (i in 1:length(keyword.list)) {
      keep[i] <- ((rowvec %*% (present.genes %in% keyword.list[[i]])) >= min.size) &
                 ((rowvec %*% (present.genes %in% keyword.list[[i]])) <= max.size)
    }
    keyword.list <- keyword.list[keep==1]
    C.names <- names(keyword.list)
    
    keyword.list <- lapply(keyword.list,function(x) x[x %in% present.genes])
    keyword.list <- lapply(keyword.list, unique)
    num.cats.vec <- sapply(keyword.list, length)
    ja <- unname(unlist(lapply(keyword.list, match, present.genes)))
    ra <- rowvec[ja]
    ia <- unname(cumsum(c(1,num.cats.vec)))
    
    dimension = c(length(C.names), length(present.genes))
    C.matrix.csr = new("matrix.csr", ra = ra, ja = as.integer(ja),
       ia = as.integer(ia), dimension = as.integer(dimension))
    C.matrix.csr = t(C.matrix.csr)
}
cat(paste("\n ",length(C.names),"catgories formed\n"))

#### 4) Prune identical categories 
  if(prune){
     drop <- apply(as.matrix(C.matrix.csr)>0,2,which)
     drop <- duplicated(sapply(drop,paste,collapse=","))
     C.matrix.csr = C.matrix.csr[,!drop]
     C.names <- C.names[!drop]
     cat(paste("  Pruned to",sum(!drop),"catgories \n"))
  }

#### 5) Return object
  if(as.matrix) {
    C.mat <- as.matrix(C.matrix.csr)
    dimnames(C.mat) <- list(present.genes,C.names)
    return(C.mat)
  } else return(list(C.mat.csr=C.matrix.csr, row.names = present.genes, col.names = C.names))
}
