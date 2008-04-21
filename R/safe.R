`safe` <-
function(X.mat, y.vec, C.mat = NULL, platform = NULL, annotate = NULL, Pi.mat = NULL, 
         local="default", global = "Wilcoxon", args.local = NULL, 
         args.global = list(one.sided=FALSE), error = "none", alpha = NA, 
         method = "permutation", min.size = 2, max.size = Inf, ...){

  if(is(X.mat,"exprSet")) {
    require("Biobase")
    pData <- pData(X.mat)
    X.mat <- exprs(X.mat)
    if (length(y.vec) == 1){
      if(is.character(y.vec)){
       if(y.vec %in% names(pData)) y.vec <- pData[[y.vec]] else {
         stop(paste("y.vec = '", y.vec, "' is not a column in pData(X.mat)",sep=""),call.=FALSE) }
      } else { 
       if(y.vec <= dim(pData)[[2]]) y.vec <- pData[[y.vec]] else {
         stop(paste("y.vec =",y.vec, " is not a column # in pData(X.mat)"),call.=FALSE) }
      }
    }
  }
  if (length(y.vec)!=dim(X.mat)[[2]]) {
    stop("Dimensions of X.mat and y.vec do not conform",call.=FALSE)}
  X.mat <- as.matrix(X.mat)

  if(is.null(C.mat)){
    if(is.null(platform) | is.null(annotate)){
      stop("C.mat or platform&annotate must be specified",call.=FALSE) 
    } else {
      require(platform,character.only=TRUE)
      names <- names(as.list(get(paste(platform,"ACCNUM",sep=""))))
      if (is.null(dimnames(X.mat)[[1]]) | (sum(dimnames(X.mat)[[1]] %in% names)==0) ) {
        stop(paste("row.names of X.mat do not conform with the '",platform,"' platform",sep=""),call.=FALSE)}
      if(substr(annotate[1],1,3)=="GO."){
        cat(paste("Building ",annotate," categories from ",platform,"GO2ALLPROBES\n",sep=""))
        C.mat <- getCmatrix(keyword.list = as.list(get(paste(platform,"GO2ALLPROBES",sep=""))), 
                            present.genes = dimnames(X.mat)[[1]], min.size = min.size,
                            max.size = max.size, GO.ont = substr(annotate,4,5))        
        C.names <- C.mat$col.names
        C.mat <- C.mat$C.mat.csr
      } else if(annotate=="KEGG"){
        cat(paste("Building ",annotate," categories from ",platform,"PATH\n",sep=""))
        C.mat <- getCmatrix(gene.list = as.list(get(paste(platform,"PATH",sep=""))), 
                            present.genes = dimnames(X.mat)[[1]], min.size = min.size,
                            max.size = max.size)
        C.names <- paste("KEGG:",C.mat$col.names,sep="")
        C.mat <- C.mat$C.mat.csr
      } else if(annotate=="PFAM"){
        cat(paste("Building ",annotate," categories from ",platform,"PFAM\n",sep=""))
        C.mat <- getCmatrix(gene.list = as.list(get(paste(platform,"PFAM",sep=""))), 
                            present.genes = dimnames(X.mat)[[1]], min.size = min.size,
                            max.size = max.size)
        C.names <- paste("PFAM:",substr(C.mat$col.names,3,100),sep="")
        C.mat <- C.mat$C.mat.csr
      } else stop(paste("Annotate = '",annotate,"' not recognized",sep=""),call.=FALSE)
   }
  } else if(class(C.mat)=="matrix"){
    require(SparseM)
    C.names <- dimnames(C.mat)[[2]]
    C.mat <- as.matrix.csr(C.mat)
    if (sum(dimnames(X.mat)[[1]]!=dimnames(C.mat)[[1]])>0) {
      cat("Warning: gene labels do not match between X.mat and C.mat")}
  } else {
    C.names <- C.mat$col.names
    C.mat <- C.mat$C.mat.csr
  }

  if (dim(C.mat)[[1]]!=dim(X.mat)[[1]]) {
    stop("Dimensions of X.mat and C.mat do not conform",call.=FALSE)}
  
  num.cats  <- dim(C.mat)[[2]]
  num.genes <- dim(X.mat)[[1]]
  num.arrays<- dim(X.mat)[[2]]

  if (local=="default") {
    if (length(unique(y.vec)) == 2){
      local <- "t.Student" 
    }  else if(class(y.vec)=="character") {
      local <- "f.ANOVA"
    } else  local <-  "t.LM"
  }
  
  if(is.null(Pi.mat)) {
    if(method %in% c("bootstrap","bootstrap.t")) Pi.mat <- 200 else Pi.mat <- 1000
  }
  if(length(unlist(Pi.mat))==1) { if(Pi.mat > 1){
    n <- length(y.vec)
    num.perms <- Pi.mat
    if(local %in% c("t.paired")){ 
      if(method!="permutation") count <- factorial(n/2) else count <- 2^(n/2-1)
      Pi.mat <- getPImatrix(block.vec = y.vec, K = Pi.mat, method = method)
    } else if(local %in% c("t.LM","f.GLM","z.COXPH")) {
      if(method!="permutation") count<-sum(choose(n,k=2:n)*choose(n-1,k=1:(n-1))) else {
         count <- factorial(n) }
      Pi.mat <- getPImatrix(n = length(y.vec), K = Pi.mat, method = method)
    } else {
      n2 <- table(y.vec)
      n3 <- 1
      for(i in 1:length(n2)) n3<-n3*sum(choose(n2[i],k=2:n2[i])*choose(n2[i]-1,k=1:(n2[i]-1)))
      if(method!="permutation") count<-n3  else count<-exp(lgamma(n+1) - sum(lgamma(n2+1)))
      Pi.mat <- getPImatrix(y.vec = y.vec, K = Pi.mat, method = method)
    }
    if(count<num.perms) cat(paste("Warning: only",round(count),"resamples exist\n"))
  } else num.perms <- 1 } else {
   num.perms <- dim(Pi.mat)[[1]]
   if(length(y.vec)!=dim(Pi.mat)[[2]]) {
     stop("Dimensions of Pi.mat and y.vec do not conform",call.=FALSE)}
  }  
 
#  print(Pi.mat[1:3,])
   
  local.stat <- get(paste("local",local,sep="."))(X.mat,y.vec,args.local) 
  u.obs <- local.stat(data = X.mat)
  names(u.obs) <- dimnames(X.mat)[[1]]

  if(!is.logical(args.global$one.sided)) stop("args.global$one.sided is missing or incorrect",call.=FALSE)
  global.stat <- get(paste("global",global,sep="."))(C.mat,u.obs,args.global)

  v.obs <- global.stat(u.obs)
  names(v.obs) <- C.names

  if(length(unlist(Pi.mat))==1) {
    return(new("SAFE", local=local, local.stat=u.obs, local.pval=as.numeric(rep(NA,num.genes)),
               global=global, global.stat=v.obs, global.pval=as.numeric(rep(NA,num.cats)),
               error=error, global.error=as.numeric(rep(NA,num.cats)), alpha=as.numeric(alpha),
               C.mat=C.mat, method=method))
  } 

  if(method== "permutation"){
    u.pvalue <- rep(1/num.perms, num.genes)
    names(u.pvalue) <- dimnames(X.mat)[[1]]

    if(error == "none"){
      emp.p <- rep(1/num.perms,num.cats)
      names(emp.p) <- C.names 
    } else {
      V.mat <- matrix(0,num.perms,num.cats)
      V.mat[1,] <- v.obs
    } 

    for(i in 2:num.perms){
      u <- local.stat(data = X.mat[,Pi.mat[i,]], resample = Pi.mat[i,])
      u.pvalue <- u.pvalue + (abs(u) >= abs(u.obs)) / num.perms

      v <- global.stat(u)
      if(error == "none") emp.p <- emp.p + (v>=v.obs)/num.perms else V.mat[i,] <- v 
      if (trunc(i/100)==i/100) cat(paste(i,"permutations completed\n"))  
    }

    if((global=="Fisher") & (!is.null(args.global$genelist.cutoff))){
        size <- (rep(1,length(u.obs)) %*% C.mat)[1,]
        if(args.global$one.sided) L <- sum(u.obs >= args.global$genelist.cutoff) else {
                                L <- sum(abs(u.obs) >= args.global$genelist.cutoff)}
        v.obs <- 1+ qhyper(v.obs,size,num.genes-size,L)
    }

    if (error != "none"){
      P.mat <- (num.perms + 1 - apply(V.mat,2,rank,ties.method="min")) / num.perms
      dimnames(P.mat)[[2]] <- dimnames(C.mat)[[2]]  

      error.p<-get(paste("error",error,sep="."))(P.mat)
      names(error.p) <- C.names
      if(is.na(alpha)) alpha <- 0.1

      return(new("SAFE",local=local, local.stat=u.obs, local.pval=u.pvalue, global=global,
             global.stat=v.obs, global.pval=P.mat[1,,drop=TRUE], error=error, alpha=alpha,
             global.error=error.p, C.mat=C.mat, method=method))
    } else {
      if(is.na(alpha)) alpha <- 0.05
      return(new("SAFE",local=local, local.stat=u.obs, local.pval=u.pvalue, global=global,
             global.stat=v.obs, global.pval=emp.p, error=error, alpha=alpha,
             global.error=as.numeric(rep(NA,num.cats)), C.mat=C.mat, method=method))
    }
 
  } else if(method=="bootstrap" | method=="bootstrap.t" |  method=="bootstrap.q"){

    if(local %in% c("f.GLM")){
      stop(paste("local = \"", local,"\" can not be used in the bootstrap"),call.=FALSE)}
    if(!global %in% c("Wilcoxon","Pearson","AveDiff")){
      stop(paste("global = \"", global,"\" cant be used in the bootstrap"),call.=FALSE)}
    if(!error %in% c("none","FWER.Bonf","FWER.Holm","FDR.BH")){
      cat(paste("WARNING: error = \"",error,"\" not available in bootstrap\n",sep=""))
      error="none"
    }
 
    if(is.null(args.local$boot.test)){
      u.pvalue <-  rep(NA, num.genes)
    } else if(args.local$boot.test=="q"){
      u.pvalue <-  rep(1/num.perms, num.genes)
    } else if(args.local$boot.test=="t"){
      u.pvalue <-  rep(NA, num.genes)
      u.sum <- u.obs 
      u2.sum <- u.obs^2
      null.local <- 0
    } else {
      cat(paste("WARNING: args.local$boot.test = \"",args.local$boot.test,
                "\" not recognized\n",sep=""))
      u.pvalue <-  rep(NA, num.genes)
      args.local$boot.test <- NULL
    }
    names(u.pvalue) <- dimnames(X.mat)[[1]]
     
    emp.p <- rep(1/num.perms,num.cats)
    names(emp.p) <- C.names

    if(global=="Wilcoxon"){
      C.size <- (rep(1,num.genes) %*% C.mat)[1,]
      null.global <- (num.genes + 1) * C.size / 2 
    } else null.global = 0

    v.sum <- v.obs 
    v2.sum <- v.obs^2

    for(i in 2:num.perms){
      u <- local.stat(data = X.mat[,Pi.mat[i,]], vector = y.vec[Pi.mat[i,]], resample = Pi.mat[i,])
      if(!is.null(args.local$boot.test)){
        if(args.local$boot.test=="q"){
          u.pvalue <- u.pvalue + (u*sign(u.obs) <= 0)/num.perms
        } else {u.sum <- u + u.sum ; u2.sum <- u^2 + u2.sum}
      }
      v <- global.stat(u)
      v.sum <- v + v.sum
      v2.sum <- v^2 + v2.sum
      emp.p <- emp.p + (v <= null.global)/num.perms
      if (trunc(i/100)==i/100) cat(paste(i,"bootstrap resamples completed\n"))
    }
    
    
    if(method=="bootstrap" | method=="bootstrap.t"){
      emp.p <-  1 - pt(((v.sum / num.perms) - null.global) / sqrt((v2.sum - v.sum^2 / 
                num.perms) / (num.perms - 1)), df = num.arrays - 1)
    } 
    
    if(!is.null(args.local$boot.test)) if(args.local$boot.test=="t") {
      u.pvalue <- 1- pt(abs(u.sum / num.perms) / sqrt((u2.sum - u.sum^2 / 
                  num.perms) / (num.perms - 1)) ,df = num.arrays - 1)
    }
    if (error == "none"){
      if(is.na(alpha)) alpha <- 0.05
      return(new("SAFE",local=local, local.stat=u.obs, local.pval=as.numeric(u.pvalue), global=global,
             global.stat=v.obs, global.pval=emp.p, error=error, alpha=alpha,
             global.error=as.numeric(rep(NA,num.cats)), C.mat=C.mat, method=method))
    } else {
      error.p<-get(paste("error",error,sep="."))(t(emp.p))
      names(error.p) <- C.names
      if(is.na(alpha)) alpha <- 0.1
      return(new("SAFE",local=local, local.stat=u.obs, local.pval=as.numeric(u.pvalue), global=global,
             global.stat=v.obs, global.pval=emp.p, error=error, alpha=alpha,
             global.error=error.p, C.mat=C.mat, method=method))
    }
  }
}

