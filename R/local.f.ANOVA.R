`local.f.ANOVA` <-
function(X.mat,y.vec,...){
  unique <- unique(y.vec)
  ANOVA.mat<-rep(1,length(y.vec))
  for(i in 2:length(unique)) {
    ANOVA.mat<-cbind(ANOVA.mat,(y.vec==unique[i])*1)
  }
  ANOVA.inner <- solve(t(ANOVA.mat) %*% ANOVA.mat)
  ANOVA.contr <- cbind(rep(0,length(unique)-1),diag(length(unique)-1))

  return(function(data = X.mat,a.x = as.matrix(ANOVA.mat),a.xpxi = ANOVA.inner,
                  a.c = as.matrix(ANOVA.contr),...){
    n<-dim(a.x)[[1]]; r<-dim(a.x)[[2]]
    a<-dim(a.c)[[1]]; m<-dim(data)[[1]]
    a.y<-t(data)
    mse<-t(rep(1,n))%*%(a.y*((diag(n)-a.x%*%a.xpxi%*%t(a.x))%*%a.y))/(n-r)
    thetahat<- a.c %*% a.xpxi %*% t(a.x) %*% a.y
    msa<-t(rep(1,a))%*%(thetahat*((solve(a.c%*%a.xpxi%*%t(a.c))%*%thetahat)))/a
    return(as.numeric(msa/mse))
  })
}

