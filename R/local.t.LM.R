"local.t.LM" <-
function(X.mat,y.vec){
  ANOVA.mat<-rep(1,length(y.vec))
  ANOVA.mat<-cbind(ANOVA.mat,y.vec)
  
  ANOVA.inner <- solve(t(ANOVA.mat) %*% ANOVA.mat)
  ANOVA.contr <- t(0:1)

  return(function(X.mat,a.x = ANOVA.mat,a.xpxi = ANOVA.inner,
                  a.c = ANOVA.contr){
    n<-dim(a.x)[[1]]; r<-dim(a.x)[[2]]
    a<-dim(a.c)[[1]]; m<-dim(X.mat)[[1]]
    a.y<-t(X.mat)
    mse<-t(rep(1,n))%*%(a.y*((diag(n)-a.x%*%a.xpxi%*%t(a.x))%*%a.y))/(n-r)
    thetahat<- a.c %*% a.xpxi %*% t(a.x) %*% a.y
    msa<-t(rep(1,a))%*%(thetahat*((solve(a.c%*%a.xpxi%*%t(a.c))%*%thetahat)))/a
    return(as.numeric(sqrt(msa/mse)*sign(thetahat)))
  })
}

