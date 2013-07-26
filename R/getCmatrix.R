getCmatrix <-
function (keyword.list = NULL, gene.list = NULL, 
          present.genes = NULL,  min.size = 2, max.size = Inf,
          by.gene = FALSE, gene.names =  NULL,
          prune = FALSE, as.matrix = FALSE,GO.ont=NULL,...)
{
##    ...: Allows arguments from versions 1.0 and 2.0 to be ignored
## keyword.list, gene.list: Objects of class 'list' containing category assigment by
##        column or by row respectively
## present.genes:      Optional character vector of row labels to match C.mat to X.mat
## GO.ont:             Optional character string, "CC","BP", or "MF" to specify ontology
## min.size, max.size: Optional numeric thresholds for gene-set size
## by.gene:            Logical to 'average' multiple probesets to a single gene should
## gene.names:         Character vector of gene names required for by.gene = TRUE
## prune:              Logical to remove redundant categories
## as.matrix:          Logical to return as class 'matrix' instead of 'SparseM'

require(SparseM)

#### 0) check for errors
  if((is.null(keyword.list) + is.null(gene.list))!=1)
    stop("One and only one of keyword.list or gene.list can be specified",call. = FALSE)

  if(is.null(keyword.list)){
    if(is.null(names(gene.list))) stop("gene.list must have names",call. = FALSE)
    if(is.null(present.genes)) present.genes <- names(gene.list)
    if(sum(present.genes %in% names(gene.list)) == 0)
      stop("No present.genes match the names of gene.list",call. = FALSE)   
  } else {
    if(is.null(present.genes)){
      cat("WARNING: rows are sorted elements of keyword.list\n")
      present.genes <- sort(unique(unlist(keyword.list)))
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
    genes.match <- match(present.genes, names(gene.list))
    present.list <- gene.list[genes.match]
    
    ignore <- sapply(present.list,is.null)
    present.list[!ignore] <- lapply(present.list[!ignore],function(x) x[!is.na(x)])
    present.list[!ignore] <- lapply(present.list[!ignore],unique)

    C.names <- sort(unique(unlist(present.list)))
    num.cats.vec <- sapply(present.list, length)

    C.matrix.csr <- new("matrix.csr",
      ra = unlist(mapply(rep, rowvec, num.cats.vec)),
      ja = match(unlist(present.list), C.names),
      ia = as.integer(cumsum(c(1,num.cats.vec))),
      dimension = as.integer(c(length(present.list),length(C.names))))
    
    size <- rep(1,length(rowvec)) %*% C.matrix.csr
    C.matrix.csr <- C.matrix.csr[, (size >= min.size) & (size <= max.size)]
    C.names <- C.names[(size >= min.size) & (size <= max.size)]
} else {
#### 3) Make C from keyword.list
    if(!is.null(GO.ont)){      
      require(GO.db)
      keep <- names(keyword.list) %in% names(as.list(GOTERM))
      keep2 <- sapply(mget(names(keyword.list)[keep], GOTERM), 
                      Ontology) %in% GO.ont
      keyword.list <- keyword.list[keep][keep2]
    }
    
    keyword.list <- lapply(keyword.list,function(x) x[x %in% present.genes])
    keyword.list <- lapply(keyword.list, unique)

    rowvec2 <- rowvec[!duplicated(present.genes)]
    present.genes2 <- present.genes[!duplicated(present.genes)]
    C.names <- names(keyword.list)
    num.cats.vec <- sapply(keyword.list, length)
    ja <- match(unlist(keyword.list), present.genes2)
    C.matrix.csr <- new("matrix.csr",
      ra = as.numeric(rowvec2)[ja],
      ja = ja,
      ia = as.integer(cumsum(c(1,num.cats.vec))),
      dimension = as.integer(c(length(C.names),length(present.genes2))))
    C.matrix.csr <- C.matrix.csr[,match(present.genes,present.genes2)]
    size <- as.numeric(C.matrix.csr %*% rep(1,length(rowvec)))
    keep <- (size >= min.size) & (size <= max.size)
    C.names <- C.names[keep]
    C.matrix.csr <- t(C.matrix.csr[keep])
    if(!is.null(col.desc)) col.desc <- col.desc[keep]
}
cat(paste(length(C.names),"categories formed\n"))

#### 4) Prune identical categories 
  if(prune){
     string <- t(C.matrix.csr)
     string <- mapply(function(s1,s2,x=string@ja) x[s2:s1],
                       string@ia[-1]-1,string@ia[-length(string@ia)])
     string <- sapply(string,paste,collapse=",")
     drop <- which(duplicated(string))

     C.matrix.csr <- C.matrix.csr[,-drop]
     C.names.drop <- C.names[drop]
     C.names <- C.names[-drop]
     if(!is.null(col.desc)) col.desc <- col.desc[-drop]
     
     col.synonym <- rep("",length(C.names))
     point <- match(string[drop],string[!drop])
     for(u in unique(point)) col.synonym[u] <-
       paste(C.names.drop[point==u],collapse="|")
     cat(paste("  Pruned to",length(C.names),"categories \n"))
  }

#### 5) Return object
  if(as.matrix) {
    C.mat <- as.matrix(C.matrix.csr)
    dimnames(C.mat) <- list(present.genes,C.names)
    return(C.mat)
  } else {
    l <- list(C.mat.csr=C.matrix.csr, row.names = present.genes, col.names = C.names)
    if(prune) l <- c(l,col.synonym=list(col.synonym))
    return(l)
  }
}
