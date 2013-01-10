sigfig <-
function(t,d){
    old.o <- options(scipen=4)                  ## Fixes default scientific notation
    v1 <- format(round(t,d),digits=d+1,trim=T)  ## Fixes digits for sig2, and any leading-white-space
    v2 <- format(signif(t,1),scientific=T)
    swap <- abs(t)<=10^(-d)
    v1[swap] <- v2[swap]
    options(scipen=old.o$scipen)
    return(v1)
}
