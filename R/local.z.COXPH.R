`local.z.COXPH` <-
function(X.mat,y.vec,args.local=NULL){

  require(survival)

  if(is.null(args.local$censor)) stop("Censoring indicators required",call.=FALSE)
  
  n = length(y.vec)

  if(is.null(args.local$covariate)){
    return(function(data = X.mat, times=y.vec, cens=args.local$censor, resample = 1:n, ...){
      m <- dim(data)[[1]]
      z <- rep(0,m)
      y <- Surv(times[resample],cens[resample])
      for(i in 1:m){
        fit <- coxph(y ~  data[i,])
        z[i] <- fit[[1]]/sqrt(fit[[2]])
      }
      return(z)
    })

  } else return(function(data = X.mat, times=y.vec, cens=args.local$censor, 
                Cov = args.local$covariate,resample = 1:n, ...){
      m <- dim(data)[[1]]
      z <- rep(0,m)
      y <- Surv(times[resample],cens[resample])
      Cov <- data.frame(Cov[resample,])
      xnam <- dimnames(Cov)[[2]]
      fmla <- as.formula(paste("y ~ data[i,] + ", paste(xnam, collapse= "+")))
      for(i in 1:m){
        fit <- coxph(fmla,data=Cov)
        z[i] <- fit[[1]][1]/sqrt(fit[[2]][1,1])
      }
      return(z)
    })
}

