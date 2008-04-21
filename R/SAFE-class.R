###### SAFE Class and Methods
##
## Note: slots changed, show updated; BB011507
##

require(methods)
require(SparseM)
setClass("SAFE",representation(local="character",local.stat="numeric",local.pval="numeric",
         global="character",global.stat="numeric",global.pval="numeric",error="character",
         global.error="numeric",C.mat="matrix.csr",alpha="numeric",method="character"))

setMethod("show","SAFE",
  function(object){
    pretty <- function(t,d) as.character(round(t,d)^(t>10^(-d))*signif(t,1)^(t<=10^(-d)))
    cat("SAFE results:\n")
    cat(paste("  Local:",object@local,"\n"))
    cat(paste("  Global:",object@global,"\n"))
    cat(paste("  Method:",object@method,"\n"))

    size <- (rep(1,length(object@local.stat)) %*% object@C.mat)[1,]

    if(object@global=="Wilcoxon"){
       stat.name <- "Mean.Rank"
       stat <- round(object@global.stat/size,1)
    } else if(object@global %in% c("Pearson","AveDiff")){
       stat.name <- "Z.stat.unstd"
       stat <- round(object@global.stat,3)
    } else if(object@global=="Fisher"){
       stat.name <- "Num.Reject"
       stat <- object@global.stat
    } else {
       stat.name <- "Global.Stat"
       stat <- object@global.stat
    }

    if (object@error=="none"){
      cat("\n")
      if(sum(is.na(object@global.pval))>0) {
        table <- data.frame(size,stat)
        dimnames(table) <- list(names(object@global.stat), c("Size",stat.name))
        print(table)
      } else {
        n <- sum(round(object@global.pval,4) <= object@alpha)
        if(n == 0) {
          cat(paste("No categories had p <=",round(object@alpha,4),"\n")) 
        } else {
        table <- data.frame(size,stat,
                       Emp.p=pretty(object@global.pval,4))[
                       order(object@global.pval),,drop=FALSE][1:n,,drop=FALSE]
        dimnames(table) <- list(names(object@global.stat[order(object@global.pval)][1:n]),
                                c("Size",stat.name,"Emp.pvalue"))
        print(table)
        }
      }
    } else {
      cat(paste("  Error:",object@error,"\n\n"))
      n = sum(object@global.error <= object@alpha)
      if(n == 0) {
        cat(paste("No categories were significant at", object@error,"<",round(object@alpha,4),"\n")) 
      } else {
        table <- data.frame(size,stat,
                       Emp.p=pretty(object@global.pval,4),
                       Adj.p=pretty(object@global.error,4))[
                       order(object@global.pval),,drop=FALSE][1:n,,drop=FALSE]
        dimnames(table) <- list(names(object@global.stat[order(object@global.pval)][1:n]),
                                c("Size",stat.name,"Emp.pvalue","Adj.pvalue"))
        print(table)
      }
    }			
  }
)

setMethod("[", "SAFE", 
  function(x, i, j,...,drop) {
    if (all(i %in% names(x@global.stat)) |
        all(i %in% 1:length(x@global.stat)) | 
        (is.logical(i) & (length(i) == length(x@global.stat)))) {
    if (all(i %in% names(x@global.stat))) i <- match(i,names(x@global.stat))
      x@global.stat <- x@global.stat[i] 
      x@global.pval <- x@global.pval[i]
      x@global.error <- x@global.error[i] 
      x@C.mat <- x@C.mat[,i]
      x@alpha <- 1
      return(x)
    } else  stop("unrecognized category(s)", call. = FALSE)
  }
)            


