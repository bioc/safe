getXYresiduals <-
function(X.mat,y.vec,Z.mat){
  if (length(unique(y.vec)) == 2){
   if (!all(sort(unique(y.vec)) == 0:1)) {
        cat(paste("Warning: y.vec is not (0,1), thus Group 1 ==",
            y.vec[1], "\n"))
        y.vec <- (y.vec == y.vec[1]) * 1
    }
    fit.y <- glm(y.vec ~ Z.mat, family = binomial(link = "logit"))
  } else fit.y <- lm(y.vec ~ Z.mat)
  y.star <- fit.y$residuals
  X.star <- X.mat
  for(i in 1:nrow(X.mat))
      X.star[i,] <- lm(X.mat[i,] ~ Z.mat)$residuals
  return(list(y.star = y.star, X.star = X.star,df.p = fit.y$rank - 1))
}
