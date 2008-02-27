`global.Fisher` <-
function(C.mat, u, args.global){
  if(is.null(args.global$genelist.cutoff)&is.null(args.global$genelist.length)) {
    stop("args.global$genelist.cutoff or args.global$genelist.length must be specified",call.=FALSE)}
  if(!is.null(args.global$genelist.cutoff) & !is.null(args.global$genelist.length)) {
    stop("Both args.global$genelist.cutoff and args.global$genelist.length cant be given",call.=FALSE)}

  m2 <-  length(u)
  size2 <- (rep(1,m2) %*% C.mat)[1,]
  if(!args.global$one.sided){
    if(!is.null(args.global$genelist.length)){
      return(function(u, C.mat2= C.mat, n  = args.global$genelist.length) {
        return(as.numeric(t(C.mat2) %*% as.numeric(rank(-abs(u))<= n)))})
    } else {
      return(function(u, C.mat2 = C.mat, cut  = args.global$genelist.cutoff, size=size2, m = m2) {
        count <- as.numeric(t(C.mat2) %*% as.numeric(abs(u) >= cut))
        return(phyper(count-1,size,m-size,sum(abs(u) >= cut)))})
    }
  } else {
    if(!is.null(args.global$genelist.length)){
      return(function(u, C.mat2 = C.mat, n  = args.global$genelist.length) {
        return(as.numeric(t(C.mat2) %*% as.numeric(rank(-u)<= n)))})
    } else {
      return(function(u, C.mat2 = C.mat, cut  = args.global$genelist.cutoff, size=size2, m = m2) {
        count <- as.numeric(t(C.mat) %*% as.numeric(abs(u) >= cut))
        return(phyper(count-1,size,m-size,sum(abs(u) >= cut)))})
    }
  }
}

