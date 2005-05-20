"safe" <-
function(X.mat, y.vec, C.mat, Pi.mat = 1000, local="default", global = "Wilcoxon",
               error = "none", write = NA, alpha = NA, method= "permutation",
               args.local = NULL, args.global = NULL){

  if(is(X.mat,"exprSet")) {
    require("Biobase")
    pData <- pData(X.mat)
    X.mat <- exprs(X.mat)
    if (length(y.vec) == 1){
      if(is.character(y.vec)){
       if(y.vec %in% names(pData)) y.vec <- pData[[y.vec]] else {
         stop("y.vec not in pData(X.mat)",call.=F) }
      } else { 
       if(y.vec <= dim(pData)[[2]]) y.vec <- pData[[y.vec]] else {
         stop("y.vec not in pData(X.mat)",call.=F) }
      }
    }
  }
  X.mat <- as.matrix(X.mat)
  C.mat <- as.matrix(C.mat)
  if (length(y.vec)!=dim(X.mat)[[2]]) {
    stop("Dimensions of X.mat and y.vec do not conform",call.=F)}
  if (dim(C.mat)[[1]]!=dim(X.mat)[[1]]) {
    stop("Dimensions of X.mat and C.mat do not conform",call.=F)}
  if (sum(dimnames(X.mat)[[1]]!=dimnames(C.mat)[[1]])>0) {
    print("WARNING: gene labels do not match between X.mat and C.mat")}
  if (length(unlist(Pi.mat))==1){
    if(Pi.mat > 1) Pi.mat <- getPImatrix(y.vec = y.vec,K = Pi.mat,method="random")
  } else  if(length(y.vec)!=dim(Pi.mat)[[2]]) {
        stop("Dimensions of Pi.mat and y.vec do not conform",call.=F)}
    
  num.cats  <- dim(C.mat)[[2]]
  if(length(unlist(Pi.mat))!=1) num.perms <- dim(Pi.mat)[[1]] else num.perms <- 1
  num.genes <- dim(X.mat)[[1]]
  num.arrays<- dim(X.mat)[[2]]
  
  if (local=="default") {
    if (length(unique(y.vec)) == 2) local <- "t.Student" else local <- "f.ANOVA" }
  if(is.null(args.local)) {
    local.stat <- get(paste("local",local,sep="."))(X.mat,y.vec) 
  } else local.stat <- get(paste("local",local,sep="."))(X.mat,y.vec,args.local) 
  u.obs <- local.stat(X.mat)
  names(u.obs) <- dimnames(X.mat)[[1]]
  u.pvalue <- rep(1/num.perms, num.genes)
  names(u.pvalue) <- dimnames(X.mat)[[1]]

  if(is.null(args.global)) {
    global.stat <- get(paste("global",global,sep="."))(C.mat,u.obs) 
  } else global.stat <- get(paste("global",global,sep="."))(C.mat,u.obs,args.global)

  v.obs <- global.stat(C.mat, u.obs)
  names(v.obs) <- dimnames(C.mat)[[2]]

  if(length(unlist(Pi.mat))==1) {
    return(new("SAFE",local=local,local.stat=u.obs,local.pval=as.numeric(rep(NA,num.genes)),global=global,
            global.stat=v.obs,global.pval=as.numeric(rep(NA,num.cats)),error=error,alpha=as.numeric(alpha),
            global.error=as.numeric(rep(NA,num.cats)),C.mat=C.mat,method=method))
  } 
  if(error == "none"){
    emp.p <- rep(1/num.perms,num.cats)
    names(emp.p) <- dimnames(C.mat)[[2]]  
  } else {
    V.mat <- matrix(0,num.perms,num.cats)
    V.mat[1,] <- v.obs
  } 
  if(!is.na(write)){
    write("Permuted Global Statistics",file=write)
    write(v.obs,file=write,append=T,ncolumns=1)
  }

  for(i in 2:num.perms){
    u <- local.stat(X.mat[,Pi.mat[i,]])
    u.pvalue <- u.pvalue + (abs(u) >= abs(u.obs)) / num.perms
    v <- global.stat(C.mat,u)
    if(error == "none") {emp.p <- emp.p + (v>=v.obs)/num.perms
      } else V.mat[i,] <- v 
    if(!is.na(write)) write(v,file=write,append=T,ncolumns=1)
    if (trunc(i/100)==i/100) print(paste(i,"permutations completed"))  
  }
  
  if (error != "none"){
    if (!is.na(write)) V.mat <- t(matrix(scan(write,sep="\n",skip=1),
                                      num.cats,num.perms))
 

    P.mat <- (num.perms + 1 - apply(V.mat,2,rank,ties.method="min")) / num.perms
    dimnames(P.mat)[[2]]<-dimnames(C.mat)[[2]]  

    error.p<-get(paste("error",error,sep="."))(P.mat)
    if(is.na(alpha)) alpha <- 0.1

    return(new("SAFE",local=local,local.stat=u.obs,local.pval=u.pvalue,global=global,
            global.stat=v.obs,global.pval=P.mat[1,],error=error,alpha=alpha,
            global.error=error.p,C.mat=C.mat,method=method))
  } else {
    if(is.na(alpha)) alpha <- 0.05
    return(new("SAFE",local=local,local.stat=u.obs,local.pval=u.pvalue,global=global,
            global.stat=v.obs,global.pval=emp.p,error=error,alpha=alpha,
            global.error=as.numeric(rep(NA,num.cats)),C.mat=C.mat,method=method))
  }

}

