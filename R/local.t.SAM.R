"local.t.SAM" <-
function(X.mat,y.vec){
  if (length(unique(y.vec))>2) {
    stop("Warning: y.vec has more than 2 elements",call.=F)
  } else if(!prod(sort(unique(y.vec))==c(0,1))) {
    print(paste("Warning: y.vec is not (0,1), thus 1 ==",y.vec[1]))
    y.vec <- (y.vec == y.vec[1])*1
  }
  SAM.stat<-function(data, vector,s0 = s0) {
    x<-as.matrix(data[,vector])
    y<-as.matrix(data[,!vector])
    n.x<-dim(x)[[2]]; n.y<-dim(y)[[2]] 
    x.m<- x %*% rep(1/n.x,n.x); y.m<- y %*% rep(1/n.y,n.y)
    ssx<-(x^2%*%rep(1,n.x)-x.m^2*n.x)
    ssy<-(y^2%*%rep(1,n.y)-y.m^2*n.y)
    t <- (x.m-y.m)/ (sqrt(ssx/n.x/(n.x-1) + ssy/n.y/(n.y-1)) + s0)
    return(t)   
  }
  n.1 <- sum(y.vec==1)
  n.0 <- sum(y.vec==0)
  mean.1 <- X.mat[,y.vec==1] %*% rep(1/n.1,n.1)
  mean.0 <- X.mat[,y.vec==0] %*% rep(1/n.0,n.0)
  ss.1 <- X.mat[,y.vec==1]^2 %*% rep(1,n.1) - mean.1^2 * n.1
  ss.0 <- X.mat[,y.vec==0]^2 %*% rep(1,n.0) - mean.0^2 * n.0
  s.vec <- sqrt((ss.1 + ss.0) * (1/n.1 + 1/n.0) / (n.1 + n.0 - 2))
  s.quant <- quantile(s.vec,probs=(0:100)/100)
  s.alpha <- rep(0, 21)
  cv.alpha <- rep(0, 21)
  for (i in 0:20){
    alpha <- i/20
    s.alpha[i+1] <- quantile(s.vec,alpha)
    v.j <- rep(0,100)
    for(j in 1:100){
      d.alpha <- SAM.stat(X.mat,y.vec==1,s0 = s.alpha[i+1])
      v.j[j] <- mad(d.alpha[(s.vec >= s.quant[j]) & (s.vec < s.quant[j+1])])
    }
    cv.alpha[i+1] <- sd(v.j) / mean(v.j)
  }
  s0.alpha <- s.alpha[order(cv.alpha)][1]
  print(paste("Alpha =", (order(cv.alpha)[1] - 1) / 20, 
              "; s_0 =", s0.alpha))
  return(function(data,vector = (y.vec == 1), s0 = s0.alpha) {
      x<-as.matrix(data[,vector])
      y<-as.matrix(data[,!vector])
      n.x<-dim(x)[[2]]; n.y<-dim(y)[[2]] 
      x.m<- x %*% rep(1/n.x,n.x); y.m<- y %*% rep(1/n.y,n.y)
      ssx<-(x^2%*%rep(1,n.x)-x.m^2*n.x)
      ssy<-(y^2%*%rep(1,n.y)-y.m^2*n.y)
      t <- (x.m-y.m)/ (sqrt(ssx/n.x/(n.x-1) + ssy/n.y/(n.y-1)) + s0)
      return(as.numeric(t))   
    })
}

