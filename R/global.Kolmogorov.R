`global.Kolmogorov` <-
function(C.mat, u,  args.global){
  m2 <-  length(u)
  size2 <- (rep(1,m2) %*% C.mat)[1,]
  if(!args.global$one.sided){
  return(function(u, C.mat2 = as.matrix(C.mat), m = m2, g.vec = size2) {
    G <- rep(1,m) %*% t(g.vec)
    ranked.Cmatrix <- C.mat2[order(-abs(u)),] * sqrt((m-G)/G) -
                      (1 - C.mat2[order(-abs(u)),]) * sqrt(G/(m-G))
    return(apply(apply(ranked.Cmatrix,2,cumsum),2,max))
  })
  } else {
  return(function(u, C.mat2 = as.matrix(C.mat),m = m2, g.vec = size2) {
    G <- rep(1,m) %*% t(g.vec)
    ranked.Cmatrix <- C.mat2[order(-u),] * sqrt((m-G)/G) -
                      (1 - C.mat2[order(-u),]) * sqrt(G/(m-G))
    return(apply(apply(ranked.Cmatrix,2,cumsum),2,max))
  })
  }
}

