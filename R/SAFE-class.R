require(methods)
setClass("SAFE",representation(local="character",local.stat="numeric",local.pval="numeric",
         global="character",global.stat="numeric",global.pval="numeric",error="character",
         global.error="numeric",C.mat="matrix",alpha="numeric",method="character"))

setMethod("show","SAFE",
  function(object){
    cat("SAFE results:\n")
    cat(paste("  Local:",object@local,"\n"))
    cat(paste("  Global:",object@global,"\n"))
    if (object@error=="none"){
      cat("\n")
      if(sum(is.na(object@global.pval))>0) {
        print(cbind(Size=apply(object@C.mat,2,sum),
              Stat=round(object@global.stat,3)))
      } else {
        n = sum(object@global.pval <= object@alpha)
        if(n == 0) {
          print(paste("No categories had p <=",round(object@alpha,4))) 
        } else {
          print(cbind(Size=apply(object@C.mat,2,sum),
                Stat=round(object@global.stat,3),
                Emp.p=round(object@global.pval,4))[
            order(object@global.pval),,drop=FALSE][1:n,,drop=FALSE])
        }
      }
    } else {
      cat(paste("  Error:",object@error,"\n\n"))
      n = sum(object@global.error <= object@alpha)
      if(n == 0) {
        print(paste("No categories were significant at", object@error,"<",round(object@alpha,4))) 
      } else {
        print(cbind(Size=apply(object@C.mat,2,sum),
              Stat=round(object@global.stat,3),
              Emp.p=round(object@global.pval,4),
              round(object@global.error,4))[order(object@global.pval),,drop=FALSE][1:n,,drop=FALSE])
      }
    }			
  }
)


setMethod("[", "SAFE", 
  function(x, i, j,...,drop) {
    if (all(i %in% names(x@global.stat)) | 
        all(i %in% 1:length(x@global.stat)) | 
        (is.logical(i) & (length(i) == length(x@global.stat)))) {
      x@global.stat <- x@global.stat[i] 
      x@global.pval <- x@global.pval[i]
      x@global.error <- x@global.error[i] 
      x@C.mat <- x@C.mat[,i,drop = FALSE]
      x@alpha <- 1
      return(x)
    } else  stop("unrecognized category(s)", call. = FALSE)
  }
)            
