`safedag` <-
function(object = NULL, ontology = NULL, top = NULL, file = NULL,
         color.cutoffs = c(0.1,0.01,0.001), filter = 0,max.GOnames = 200){
  require("GOstats");require("Rgraphviz")

  if(is.null(ontology)) {
    ontology <- Ontology(mget(names(object@global.stat)[1],GOTERM)[[1]])
    cat(paste("Ontology unspecified so assumed to be GO.",ontology,"\n",sep=""))
  } else if(!ontology %in% c("GO.CC","GO.BP","GO.MF")) {
    stop("ontology must be \"GO.CC\", \"GO.BP\", or \"GO.MF\"", call.=FALSE)
    ontology = substr(ontology,4,5)
  }

  if(!is.null(top)){
    keep2 <- top
    stop <- 0
    names <- names(as.list(get(paste("GO",ontology,"CHILDREN",sep=""))))
    while(stop==0){
      parents <- keep2[keep2 %in% names]
      children <- unlist(mget(parents,get(paste("GO",ontology,"CHILDREN",sep=""))))
      stop <- all(children %in% keep2)
      keep2 <- unique(c(keep2,children))
    }
  } else {
    top <- c("GO:0005575","GO:0008150","GO:0003674")[c("CC","BP","MF") %in% ontology]
    keep2 <- NULL
  }
  parents <- get(paste("GO",ontology,"PARENTS",sep=""))

  C.names <- names(object@global.stat)
  keep <- C.names %in% names(as.list(parents)) 

  cat(paste(sum(keep)," of ",length(C.names)," categories in GO.",
                     ontology,"\n",sep=""))

  if(!is.null(keep2)){
    keep[!C.names %in% keep2] <- FALSE
    cat(paste(sum(keep),"of",length(C.names),"categories below",top,"\n"))
  }

  p <- object@global.pval
  cat(paste("\n",sum(p<=color.cutoffs[3]),"categories with p <=",color.cutoffs[3],"\n"))
  if(filter<3) cat(paste("",sum(p<=color.cutoffs[2]),"categories with p <=",color.cutoffs[2],"\n"))
  if(filter<2) cat(paste("",sum(p<=color.cutoffs[1]),"categories with p <=",color.cutoffs[1],"\n"))

  if(filter) keep[p>color.cutoffs[filter]] <- FALSE

  color <- 8 - (p<=color.cutoffs[1])*6 +
               (p<=color.cutoffs[2]) + 
               (p<=color.cutoffs[3])

  
  g <- GOGraph(C.names[keep], parents)
  nda <- makeNodeAttrs(g)
  nodes <- names(nda$fillcolor)
  if(is.null(keep2)) g <- subGraph(nodes[!nodes %in% c("all","top")],g) else
                     g <- subGraph(nodes[nodes %in% keep2],g) 
  nda <- makeNodeAttrs(g)
  nodes <- names(nda$fillcolor)
  match <- match(nodes,C.names[keep], nomatch=0)

  nda$fillcolor[match] <- color[keep]
  nda$fillcolor[!nodes %in% C.names] <- "white"
  if(length(nodes)>max.GOnames) nda$label[nda$label!=""] <- " " else { 
    if(is.null(keep2)) nda$label[nodes==top] <- ontology }
  g <- agopen(g, nodeAttrs = nda,name="whatever")
   
  if(!is.null(file)) postscript(file) #else x11()
  par(mfrow=c(1,1))
  plot(g)
  if(!is.null(file)) dev.off()
}

