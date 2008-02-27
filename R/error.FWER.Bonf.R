`error.FWER.Bonf` <-
function(P.matrix){
  cap <- (P.matrix[1,]*dim(P.matrix)[[2]]) <= 1
  return((P.matrix[1,]*dim(P.matrix)[[2]])^cap)
}

