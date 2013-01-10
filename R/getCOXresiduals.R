getCOXresiduals <-
function(X.mat,time,cens,Z.mat=NULL,method="martingale"){
  require(survival)
  if(is.null(Z.mat)) Z.mat <- cbind(rep(1,ncol(X.mat)))
  surv <- Surv(time,cens)
  fit.z <- coxph(surv ~ Z.mat)
  y.star <- residuals(fit.z, type=method)
  X.star <- X.mat
  for(i in 1:nrow(X.mat))
      X.star[i,] <- lm(X.mat[i,] ~ Z.mat)$residuals
  return(list(y.star = y.star, X.star = X.star))
}
