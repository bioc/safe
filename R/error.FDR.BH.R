`error.FDR.BH` <-
function(P.matrix){
  num.cat<-dim(P.matrix)[[2]]
  order <- order(order(P.matrix[1,]))
  BH <- P.matrix[1,order(P.matrix[1,])] * 
        num.cat/(1:num.cat)
  for(i in (num.cat-1):1) BH[i]<-min(BH[i],BH[i+1])
  return(BH[order])
}

