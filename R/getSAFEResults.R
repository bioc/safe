`getSAFEResults` <- 
function(safe.obj, alpha=NULL){
  if(is.null(alpha)) alpha <- safe.obj@alpha 
  if (safe.obj@error=="none"){  
    keep <- safe.obj@global.pval <= alpha
  } else keep <- safe.obj@global.error <= alpha

  if(safe.obj@global=="Wilcoxon"){
    global <- "Wilcoxon Rank Sum"
  } else if(safe.obj@global=="AveDiff"){
    global <- "Average Difference (Welch t)"
  } else if(safe.obj@global=="Pearson"){
    global <- "Pearson Chi-square"
  } else if(safe.obj@global=="Fisher"){
    global <- "Fishers Exact P-value"
  } else global <- safe.obj@global
 
  if(safe.obj@local=="t.Student"){
    local <- "Student's t-statistic"
  } else if(safe.obj@local=="t.Welch"){
     local <- "Welch t-statistic"
  } else if(safe.obj@local=="t.paired"){
     local <- "Paired t-statistic"
  } else if(safe.obj@local=="z.COXPH"){
     local <- "Cox PH z-statistic (Wald)"
  } else if(safe.obj@local=="f.ANOVA"){
     local <- "ANOVA F-statistic"
  }  else if(safe.obj@local=="t.LM"){
     local <- "Linear regression t-statistic"
  }  else local <- safe.obj@local

  if(safe.obj@error=="FDR.BH"){
    error <- "FDR (Benjamini-Hochberg)"
  } else if(safe.obj@error=="FDR.YB"){
    error <- "FDR (Yekutieli-Benjamini)"
  } else if(safe.obj@error=="FWER.Holm"){
    error <- "FWER (Holm's)"
  } else if(safe.obj@error=="FWER.Bonf"){
    error <- "FWER (Bonferroni)"
  } else if(safe.obj@error=="FWER.WY"){
    error <- "FWER (Westfall-Young)"
  } else error <- safe.obj@error

  size <- (rep(1,length(object@local.stat)) %*% object@C.mat)[1,]

  return(list(local = local,                            # character-string describing local statistic
              gene.names = names(safe.obj@local.stat),  # character-vector of gene-names
              gene.stats = safe.obj@local.stat,         # numeric vector of gene-specific statistics
              gene.pvals = safe.obj@local.pval,         # numeric vector of gene-specific p-values
              global = global,                          # character-string describing global statistic
              cat.names = names(safe.obj@global.stat)[keep],            # character-vector of cat-names
              cat.sizes = size[keep],                         # Numeric vector of category sizes 
              cat.stats = safe.obj@global.stat[keep],         # numeric vector of cat-specific statistics
              cat.pvals = safe.obj@global.pval[keep],         # numeric vector of cat-specific p-values
              cat.error = safe.obj@global.error[keep],        # numeric vector of cat-specific adjusted p-values
              error = error,                                  # character-string describing error adj. method
              method = safe.obj@method,                       # character-string describing resampling method
              alpha = alpha                              # numeric value giving cutoff for cat-significance
              )
          )
}  
