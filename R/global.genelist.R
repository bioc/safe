"global.genelist" <-
function(C.mat,u, args.global){
  one.sided <- args.global[[1]]
  cut.off <- args.global[[2]]
  if(one.sided) {
    num.genelist <- sum(u>=cut.off)
    return(function(C.mat,u, n  = num.genelist) {
      return(as.numeric(t(C.mat) %*% (rank(-u)<=n)))
    })
  } else {
    num.genelist <- sum(abs(u)>=cut.off)
    return(function(C.mat,u, n  = num.genelist) {
      return(as.numeric(t(C.mat) %*% (rank(-abs(u))<=n)))
    })
  }
}

