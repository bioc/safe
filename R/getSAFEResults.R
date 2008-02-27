`getSAFEResults` <- 
function(safe.obj, alpha=NULL){
  if(is.null(alpha)) alpha <- safe.obj@alpha 
  if (safe.obj@error=="none"){  
    keep <- safe.obj@global.pval <= alpha
  } else keep <- safe.obj@global.error <= alpha
  return(list(local = safe.obj@local,                   # character-string describing local statistic
              gene.names = names(safe.obj@local.stat),  # character-vector of gene-names
              gene.stats = safe.obj@local.stat,         # numeric vector of gene-specific statistics
              gene.pvals = safe.obj@local.pval,         # numeric vector of gene-specific p-values
              global = safe.obj@global,                 # character-string describing global statistic
              cat.names = names(safe.obj@global.stat)[keep],  # character-vector of cat-names
              cat.stats = safe.obj@global.stat[keep],         # numeric vector of cat-specific statistics
              cat.pvals = safe.obj@global.pval[keep],         # numeric vector of cat-specific p-values
              cat.error = safe.obj@global.error[keep],        # numeric vector of cat-specific adjusted p-values
              error = safe.obj@error,                   # character-string describing error adj. method
              method = safe.obj@method,                 # character-string describing resampling method
              alpha = alpha                   # numeric value giving cutoff for cat-significance
              )
          )
}  
