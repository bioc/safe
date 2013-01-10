global.Fisher <-
function(C.mat, u, args.global){
  if(!is.null(args.global$genelist.cutoff)) {
    stop("args.global$genelist.cutoff is nolonger used for global = Fisher",call.=FALSE)}
  if(is.null(args.global$genelist.length)) {
    stop("args.global$genelist.length must be specified",call.=FALSE)}

  m2 <-  length(u)
  size2 <- (rep(1,m2) %*% C.mat)[1,]
  if(!args.global$one.sided){
      return(function(u, C.mat2 = C.mat, n  = args.global$genelist.length) {
        return(as.numeric(t(C.mat2) %*% as.numeric(rank(-abs(u))<= n)))})
  } else {
      return(function(u, C.mat2 = C.mat, n  = args.global$genelist.length) {
        return(as.numeric(t(C.mat2) %*% as.numeric(rank(-u)<= n)))})
  }
}
