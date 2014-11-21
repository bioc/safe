safe.toptable <-
function(safe, number=10, pretty=TRUE, description=TRUE){
##  safe: Objects of class 'SAFE'
##  number: number of top categories to return
##  pretty: whether or not to call sigfigs() for rounding
##  description: whether or not to show Description column
  require(SparseM)
  if(description){
    names <- names(safe@global.stat)
    desc <- as.character(rep(NA,length(safe@global.stat)))
    doGO <- substr(names,1,3) == "GO:"
    if(sum(doGO)){
      require(GO.db)
      desc[doGO] <- sapply(mget(names[doGO],GOTERM),Term)
    }
#    doKEGG <- substr(names,1,5) == "KEGG:"
#    if(sum(doKEGG)){
#      require(KEGG.db)
#      terms <- substr(names[doKEGG],6,10)
#      desc[doKEGG] <- unlist(mget(terms,KEGGPATHID2NAME))
#    }
    doPFAM <- substr(names,1,5) == "PFAM:"
    if(sum(doPFAM)){
       require(PFAM.db)
       terms <- gsub("PFAM:","PF",names[doPFAM])
       p.lists <- mget(terms,PFAMSCOP,ifnotfound=NA)
       flatten <- function(c.list) paste(unlist(c.list), collapse = ", ")
       desc[doPFAM] <- unlist(lapply(p.lists, flatten))
    }
    doREACT <- substr(names,1,8) == "REACTOME"
    if(sum(doREACT)){
      require(reactome.db)
      terms <- gsub("REACTOME:","",names[doREACT])
      ## use mapIds() instead of mget()
      tmpRes <- mapIds(reactome.db, keys=terms, column='PATHNAME',
                       keytype='PATHID', multiVals='list')
      desc[doREACT] <- sapply(tmpRes,
                          function(x) strsplit(x[1],": ")[[1]][2])
      ## desc[doREACT] <- sapply(mget(terms,reactomePATHID2NAME),
      ##                     function(x) strsplit(x[1],": ")[[1]][2])
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
