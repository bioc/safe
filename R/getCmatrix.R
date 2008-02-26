"getCmatrix" <-
function(list = NULL,vector = NULL, file = NULL, delimiter = ",", 
                     present.genes = NULL, min.size = 0, max.size = Inf, ...){
  if(is.null(list)) {  
    if(is.null(vector)) {
      vector<-as.character(read.table(file,...)[,1])
    }
    list<-strsplit(vector,delimiter)
  }
  unique<-sort(unique(unlist(list)))
  C.matrix<-matrix(FALSE,length(list),length(unique))
  dimnames(C.matrix)<-list(names(list),unique)
  for(i in 1: length(list)){
    m<-length(list[[i]])
    if (m) for(j in 1:m){
      if(!is.na(list[[i]][j])) {
        C.matrix[i,]<-C.matrix[i,]|(list[[i]][j]==unique)
      }
    }
  }
  if (!is.null(present.genes)) {
    if (!is.character(present.genes)) {
      if(!is.null(names(present.genes))) {
        if(sum(names(present.genes)!=dimnames(C.matrix)[[1]])>0) {
          stop("Names of C matrix and present.genes do not conform",call.=FALSE)
        } 
      } 
      C.matrix <- C.matrix[present.genes,]
    } else {
      C.matrix <- C.matrix[match(present.genes,dimnames(C.matrix)[[1]],nomatch=0),]
    }
  }
  size <- apply(C.matrix,2,sum)
  C.matrix<-C.matrix[,(size>=min.size)&(size <= max.size)]
  return(C.matrix*1)
}

