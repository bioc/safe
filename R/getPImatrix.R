"getPImatrix" <-
function(n = NULL, y.vec = NULL, block.vec = NULL , K, method="random"){
  if (!is.null(y.vec)) n <- length(y.vec)
  if (!is.null(block.vec)) n <- length(block.vec)
  if (length(n)>1) {
    n <- length(n)
    print(paste("Warning: n is given as a vector, and is assumed to be",n))
  }
  Pi<-matrix(NA,K,n)
  Pi[1,]<-1:n
  if(is.null(block.vec)){
    for (i in 2:K) Pi[i,]<-sample(1:n)
  } else {
   unique<-unique(block.vec)
   for (i in 2:K){
     for(j in 1:length(unique)){
       Pi[i,block.vec==unique[j]]<-sample((1:n)[block.vec==unique[j]])
     }
   }
 }
 return(Pi)
}

