`getPImatrix` <-
function(n = NULL, y.vec = NULL, block.vec = NULL , K, 
         method="permutation"){
  if (!is.null(y.vec)) n <- length(y.vec)
  if (!is.null(block.vec)) n <- length(block.vec)
  if (length(n)>1) {
    n <- length(n)
    cat(paste("Warning: n is given as a vector, and is assumed to be",n,"\n"))
  }
  Pi<-matrix(NA,K,n)
  Pi[1,]<-1:n
  
  if(substr(method,1,9)=="bootstrap"){
    if(is.null(y.vec)&is.null(block.vec)) {
       for(i in 2:K) Pi[i,]<-sample(1:n,replace=TRUE) 
    } else if(is.null(block.vec)){
      unique<-unique(y.vec)
      for (i in 2:K){
        for(j in 1:length(unique)){
          Pi[i,y.vec==unique[j]]<-sample((1:n)[y.vec==unique[j]],replace=TRUE)
        }
      }
    } else {
      pairs <- unique(abs(block.vec))
      for (i in 2:K){
        for(j in 1:length(pairs)){
          sample <- sample(pairs,1)       
          Pi[i,match(c(pairs[j],-pairs[j]),block.vec)] <- (1:n)[match(
             c(sample,-sample),block.vec)]
        }
      }
    }
    return(Pi)
  } else {
    if(is.null(block.vec)) for (i in 2:K) Pi[i,]<-sample(1:n) else {
      if(is.numeric(block.vec)) block.vec <- abs(block.vec)
      unique<-unique(block.vec)
      for (i in 2:K){
        for(j in 1:length(unique)){
          Pi[i,block.vec==unique[j]]<-sample((1:n)[block.vec==unique[j]])
        }
      }
    }
    return(Pi)
  }
}

