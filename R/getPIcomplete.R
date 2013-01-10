getPIcomplete <-
function(y.vec){
## y.vec: Object of class 'numeric', 'integer' or 'character' with binary response.
  if (length(unique(y.vec)) == 2){
   if (prod(sort(unique(y.vec)) != 0:1)) y.vec <- (y.vec == y.vec[1]) * 1
  } else stop("Exhaustive permutation is only constructed for the two-sample problem",call. = FALSE)
  rank <- order(order(y.vec))
  n <- length(y.vec); n1 <- sum(y.vec==min(y.vec))
  M <- matrix(1:n1,1,n1)
  n2 <- n1 + 1
  for(n2 in ((n1+1):n)){
    M2 <- M; M2[,1] <- n2
    for(i in 2:n1){
      temp <- M
      temp[,i] <- n2
      M2 <- rbind(M2,temp)
    }
    M <- rbind(M,t(apply(M2,1,sort)))
    drop <- duplicated(apply(M,1,paste,collapse=" "))
    M <- M[!drop,]
  }
  M <- t(apply(M,1,function(x,vec) c(x,vec[!vec %in% x]),1:n))
  M <- t(apply(M,1,function(x,r) x[r],rank))
  return(M)
}
