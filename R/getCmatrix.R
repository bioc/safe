`getCmatrix` <-
function(keyword.list = NULL, gene.list = NULL, vector = NULL, file = NULL, 
         delimiter = ",", present.genes = NULL, GO.ont = NULL, 
         min.size = 0, max.size = Inf, as.matrix = FALSE,...){

  require(SparseM)

  if(is.null(keyword.list)){
    if(is.null(gene.list)) {  
      if(is.null(vector)) {
        vector<-as.character(read.table(file,...)[,1])
      }
      gene.list<-strsplit(vector,delimiter)
    }
    C.names <- sort(unique(unlist(gene.list)))
    row.names <- names(gene.list)

    if (!is.null(present.genes)) {
      if(prod(present.genes %in% row.names) == 0) {
        cat("WARNING: Some present.genes do not appear in gene.list\n")}
      if(sum(present.genes %in% row.names) == 0) {
         stop("No present.genes matched the names of gene.list",call. = FALSE)}
    } else present.genes <- row.names  

    genes.match = match(present.genes, row.names)
    present.list <- gene.list[genes.match]
    names(present.list) = present.genes

# ra is the list of values in the csr matrix (all 1's in our case)
# ja is the column position in the matrix where each value goes for a given row.
# ia is the position in ra and ja where the next row starts.  ia has a size of
# number of rows + 1
    ra = 1
    ja = 1
    ia = 1
    raja.index = 1
    cat("  Categories completed:")
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
          ra[raja.index:last.index] = rep(1, num.cats)
          ja[raja.index:last.index] = cat.matches
          ia[g] = raja.index
          raja.index = raja.index + num.cats
       } else ia[g] = raja.index 
       } else ia[g] = raja.index
    }
    ia[g+1] = raja.index
    dimension = c(num, length(C.names))
    C.matrix.csr = new("matrix.csr", ra = as.integer(ra), ja = as.integer(ja),
       ia = as.integer(ia), dimension = as.integer(dimension))

  } else {
    if(!is.null(GO.ont)){
      require(GO.db)
      keep <- sapply(mget(names(keyword.list),GOTERM),Ontology) %in% GO.ont
      keyword.list <- keyword.list[keep]
    } 
    if(is.null(present.genes)) present.genes <- sort(unique(unlist(keyword.list)))
    for (i in 1:length(keyword.list)) {
      keep[i] <- (sum(present.genes %in% keyword.list[[i]]) >= min.size) & 
                 (sum(present.genes %in% keyword.list[[i]]) <= max.size)
    }
    keyword.list <- keyword.list[keep]
    C.names <- names(keyword.list)

# ra is the list of values in the csr matrix (all 1's in our case)
# ja is the column position in the matrix where each value goes for a given row.
# ia is the position in ra and ja where the next row starts.  ia has a size of
# number of rows + 1
    ra = 1
    ja = 1
    ia = 1
    raja.index = 1
    cat("  Categories completed:")
    num <- length(keyword.list)
    for(g in 1:num){
       if (g %in% floor(seq(num/5,num,num/5))) cat(paste(ceiling(g/num*100),"% ",sep=""))
       gene.cats = unique(keyword.list[[g]][keyword.list[[g]] %in% present.genes])
       num.cats = length(gene.cats)
       if(!is.null(gene.cats)){ if(!is.na(gene.cats[1])){
          cat.matches = match(gene.cats, present.genes)
          last.index = raja.index + num.cats - 1
          ra[raja.index:last.index] = rep(1, num.cats)
          ja[raja.index:last.index] = cat.matches
          ia[g] = raja.index
          raja.index = raja.index + num.cats
       } else ia[g] = raja.index } else ia[g] = raja.index
    }
    ia[g+1] = raja.index
    dimension = c(length(C.names), length(present.genes))
    C.matrix.csr = new("matrix.csr", ra = as.integer(ra), ja = as.integer(ja),
       ia = as.integer(ia), dimension = as.integer(dimension))
    C.matrix.csr = t(C.matrix.csr)
  }
# The SparseM matrix does not support the apply() function. So multiply
# by a row vector of 1's to get the sum of genes in each category.
  rowvec = as.matrix.csr(matrix(rep(1, dim(C.matrix.csr)[1]), nrow = 1))
  size = as.matrix(rowvec %*% C.matrix.csr)
  C.matrix.csr <- C.matrix.csr[, (size >= min.size) & (size <= max.size)]
  C.names = C.names[(size >= min.size) & (size <= max.size)]
  cat(paste("\n ",length(C.names),"catgories formed\n"))
  if(as.matrix) {
    C.mat <- as.matrix(C.matrix.csr)
    dimnames(C.mat) <- list(present.genes,C.names)
    return(C.mat)
  } else return(list(C.mat.csr=C.matrix.csr, row.names = present.genes, col.names = C.names))
}
