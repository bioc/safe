"global.Kolmogorov" <-
function(C.mat,u){
  size <- apply(C.mat,2,sum)
  return(function(C.mat,u, m = dim(C.mat)[[1]], g.vec = size) {
    G <- rep(1,m) %*% t(g.vec)
    ranked.Cmatrix <- C.mat[order(-abs(u)),] * sqrt((m-G)/G) -
                      (1 - C.mat[order(-abs(u)),]) * sqrt(G/(m-G))
    return(apply(apply(ranked.Cmatrix,2,cumsum),2,max))
  })
}

