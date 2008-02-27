`global.Pearson` <-
function(C.mat,u, args.global){
  if(is.null(args.global$genelist.cutoff)) {
    stop("args.global$genelist.cutoff must be specified",call.=FALSE)}
  m2 <-  length(u)
  size2 <- (rep(1,m2) %*% C.mat)[1,]
  if(!args.global$one.sided){
      return(function(u, C.mat2 = C.mat,cut  = args.global$genelist.cutoff, size=size2, m = m2) {
        o1 <- (t(C.mat2) %*% as.numeric(abs(u) >= cut))
        o2 <- sum(abs(u) >= cut) - o1
        diff <- (o1/size - o2/(m-size))/sqrt((o1+o2)*(1-o1/m-o2/m)/size/(m-size))
      return(as.numeric(diff))})
  } else {
      return(function(u, C.mat2 = C.mat,cut  = args.global$genelist.cutoff, size=size2, m = m2) {
        o1 <- (t(C.mat2) %*% as.numeric(u >= cut))
        o2 <- sum(u >= cut) - o1
        diff <- (o1/size - o2/(m-size))/sqrt((o1+o2)*(1-o1/m-o2/m)/size/(m-size))
      return(as.numeric(diff))})
  } 
}

