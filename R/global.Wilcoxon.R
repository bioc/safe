"global.Wilcoxon" <-
function(C.mat,u){
  return(function(C.mat,u) return(as.numeric(t(C.mat) %*% rank(abs(u)))))
}

