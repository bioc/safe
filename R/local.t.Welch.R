`local.t.Welch` <-
function(X.mat,y.vec,...){
    if (length(unique(y.vec))>2) {
       stop("local = \"t.Welch\", but y.vec has more than 2 elements",call.=FALSE)
    } else if(sum(sort(unique(y.vec))==c(0,1))!=2) {
      cat(paste("Warning: y.vec is not (0,1), thus Group 1 ==",y.vec[1],"\n"))
      y.vec <- (y.vec == y.vec[1])*1
    }
    return(function(data,vector = y.vec,...) {
      x<-as.matrix(data[,vector==1])
      y<-as.matrix(data[,!vector==1])
      n.x<-dim(x)[[2]]; n.y<-dim(y)[[2]] 
      x.m<- x %*% rep(1/n.x,n.x); y.m<- y %*% rep(1/n.y,n.y)
      ssx<-(x^2%*%rep(1,n.x)-x.m^2*n.x)
      ssy<-(y^2%*%rep(1,n.y)-y.m^2*n.y)
      t <- (x.m-y.m)/ sqrt(ssx/n.x/(n.x-1) + ssy/n.y/(n.y-1))
      return(as.numeric(t))   
    })
}

