safe.toptable <-
function(safe, number=10, annotate = c("GO","KEGG","PFAM"),
                            pretty=TRUE, description=TRUE){
##  safe: Objects of class 'SAFE'
##  number: number of top categories to return
##  annotate: Charcter vector of pathway annotation(s) for "GO" and "KEGG" categories
##  pretty: whether or not to call sigfigs() for rounding
##  description: whether or not to show Description column
  require(SparseM)
  if(description){
    names <- names(safe@global.stat)
    desc <- as.character(rep(NA,length(safe@global.stat)))
    if("GO" %in% annotate){
       require(GO.db)
       keep <- substr(names,1,3) == "GO:"
       if(sum(keep)) desc[keep] <- sapply(mget(names[keep],GOTERM),Term)
    }
    if("KEGG" %in% annotate){
       require(KEGG.db)
       keep <- substr(names,1,5) == "KEGG:"
       if(sum(keep)){
         term <- substr(names[keep],6,10)
         desc[keep] <- unlist(mget(term,KEGGPATHID2NAME))
       }
    }
    if("PFAM" %in% annotate){
       require(PFAM.db)
       keep <- substr(names,1,2) == "PF"
       if(sum(keep)){
         p.lists <- mget(names[keep],PFAMSCOP,ifnotfound=NA)
         flatten <- function(c.list) paste(unlist(c.list), collapse = ", ")
         desc <- unlist(lapply(p.lists, flatten))
       }
    }
    table <- data.frame(GenesetID = names(safe@global.stat),
                        Size = round((rep(1,length(safe@local.stat)) %*%
                                      safe@C.mat)[1,],2),
                        Statistic = safe@global.stat,
                        P.value = safe@global.pval,
                        Adj.p.value = safe@global.error,
                        Description = desc)
  } else{
    table <- data.frame(GenesetID = names(safe@global.stat),
                        Size = round((rep(1,length(safe@local.stat)) %*%
                                      safe@C.mat)[1,],2),
                        Statistic = safe@global.stat,
                        P.value = safe@global.pval,
                        Adj.p.value = safe@global.error)
  }
    table <- table[order(table$P.value,partial = 1/abs(table$Statistic)),]
    if(pretty){
      table$P.value <- sigfig(table$P.value,4)
      table$Adj.p.value <- sigfig(table$Adj.p.value,4)
    }
    rownames(table) <- 1:length(safe@global.stat)
    return(table[1:number,])
}
