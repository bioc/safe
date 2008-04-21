`global.AveDiff` <-
function(C.mat, u, args.global ){
  m <-  length(u)
  size <- (rep(1,m) %*% C.mat)[1,]
  if(!args.global$one.sided){
    return(function(u, C.mat2 = C.mat, sizeB = size, mB = m){
      m1 <- (t(C.mat2) %*% abs(u))/sizeB 
      m2 <- (sum(abs(u)) - m1*sizeB)/(mB-sizeB)
      s1 <- (t(C.mat2) %*% abs(u)^2) - m1^2 * sizeB
      s2 <- (sum(abs(u)^2) - s1 - m1^2 * sizeB) - m2^2 * (mB-sizeB)
      diff <- (m1 - m2)/sqrt(s1/sizeB/(sizeB-1) + s2/(mB-sizeB)/(mB-sizeB-1))
      return(as.numeric(diff))})
  } else return(function(u, C.mat2 = C.mat, sizeB = size ,mB = m){
      m1 <- (t(C.mat2) %*% u)/sizeB 
      m2 <- (sum(u) - m1*sizeB)/(mB-sizeB)
      s1 <- (t(C.mat2) %*% u^2) - m1^2 * sizeB
      s2 <- (sum(u^2) - s1 - m1^2 * sizeB) - m2^2 * (mB-sizeB)
      diff <- (m1 - m2)/sqrt(s1/sizeB/(sizeB-1) + s2/(mB-sizeB)/(mB-sizeB-1))
      return(as.numeric(diff))})
}

