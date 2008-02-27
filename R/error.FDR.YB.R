`error.FDR.YB` <-
function(P.matrix){
  order <- order(order(P.matrix[1,]))
  P.matrix <- P.matrix[,order(P.matrix[1,])]
  num.perms<-dim(P.matrix)[[1]]
  num.cat<-dim(P.matrix)[[2]]

  s.hat<-rep(NA,num.cat)
  YB<-rep(NA,num.cat)
  YB.funct<-function(k.vec,p,s){
    v.k<-sum(k.vec<=p)
    if ((v.k+s)>0) return(v.k/(v.k+s)) else return(0)
  }
  for(i in 1:num.cat){
    s.hat[i]<-max(0,sum(P.matrix[1,]<=P.matrix[1,i])-
                  sum(P.matrix[-1,]<=P.matrix[1,i])/(num.perms-1))
    YB[i]<-mean(apply(P.matrix[-1,],1,YB.funct,P.matrix[1,i],s.hat[i]))
  }
  for(i in (num.cat-1):1) YB[i]<-min(YB[i],YB[i+1])
  return(YB[order])
}

