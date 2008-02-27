`global.Wilcoxon` <-
function(C.mat, u, args.global){
  if(! args.global$one.sided){
    return(function(u,C.mat2=C.mat) return(as.numeric(t(C.mat2) %*% rank(abs(u)))))
  } else return(function(u,C.mat2=C.mat) return(as.numeric(t(C.mat2) %*% rank(u))))
}

