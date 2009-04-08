`safeplot` <-
function(safe = NULL, cat.name = NULL, c.vec = NULL, local.stats = NULL, p.val = NULL, 
                     one.sided = NA, limits = c(-Inf,Inf), extreme = NA, italic = FALSE , 
                     x.label = "Ranked local statistic"){
    if(!is.null(safe)){
      local.stats <- safe@local.stat
      C.names <- names(safe@global.stat)
      if (is.na(one.sided)) if(substr(safe@local,1,1)=="f") one.sided<-TRUE else one.sided<-FALSE
      if(is.null(cat.name)) cat.name <- C.names[order(safe@global.pval)][1]
      if(one.sided & length(limits)==1) limits <- c(-Inf,limits)
      if (prod(limits == c(-Inf,Inf)) == 1 & prod(1-is.na(safe@local.pval))) {
        if(!one.sided) limits[1] <- min(local.stats[(safe@local.pval > 0.05)&(local.stats < 0)])
        limits[2] <- max(local.stats[(safe@local.pval > 0.05)&(local.stats > 0)])
        cat(paste("Shaded limits are the largest local statistics with p > 0.05\n",
                   "  Limits = (",round(limits[1],3),",",round(limits[2],3),")\n"))
      }
      c.vec <- as.matrix(safe@C.mat[,C.names == cat.name])[,1]
      if(is.null(p.val)) p.val <- safe@global.pval[C.names == cat.name]
    } 
    if(!is.null(names(local.stats))) gene.names <- names(local.stats) else gene.names <- NULL
    if(is.na(extreme)) extreme <- (sum(c.vec) > 24)

    g <- sum(c.vec==1)
    m <- length(c.vec)
  
    local.sorted <- sort(local.stats)
    c.vec.sorted <- c.vec[order(local.stats)]
    if(!is.null(gene.names)) gene.names.sorted <- gene.names[order(local.stats)]
    cdf <- cumsum(c.vec.sorted)/ g

    nf <- layout(matrix(c(1,2),2,1,byrow=TRUE), c(4,4), c(1,3), TRUE)
    par(plt=c(.1,.9,0,1),adj=0)
    plot(1:m,rep(0.5,m),ylim=c(0.04,1),xaxs="i",axes=FALSE,ylab="",xlab= "",col=0)
    for(i in 1:m) if(c.vec.sorted[i]==1) points(rep(i,2),c(0.005,0.1),type="l",lwd=2)
    if(!is.null(gene.names)){
      zero <- sum(local.stats<0) + 0.5
      if (extreme) c.vec.sorted[(local.sorted>=limits[1])&(local.sorted<=limits[2])] <- 0
      tick<-(1:m)[c.vec.sorted==1]
      gA <- sum(tick<zero)
      gB <- sum(tick>zero)
      tick.adj<-c(rep(max(tick[1],0.01*m),gA),rep(min(tick[gA+gB],0.99*m),gB))
      i <- gA+gB;  stopB <- 0
      if (gB) lines(c(tick.adj[i],tick[i]), c(0.25, 0.1))
      if (gB>1) while(stopB == 0){
        i <- i - 1
        tick.adj[i]<-min(tick.adj[i+1]-0.02*m,tick[i])
        if(tick.adj[i]>m/2) lines(c(tick.adj[i],tick[i]), c(0.25, 0.1)) else {
          stopB <- i+1}
        if(i == gA+1) stopB <- gA+1
      } else stopB <- 0
      i <- 1;  stopA <- 0
      if (gA) lines(c(tick.adj[i],tick[i]), c(0.25, 0.1)) 
      if (gA>1) while(stopA == 0){
        i <- i + 1
        tick.adj[i]<- max(tick.adj[i-1]+0.02*m,tick[i])
        if(tick.adj[i]<m/2) lines(c(tick.adj[i],tick[i]), c(0.25, 0.1)) else {
          stopA <- i-1}
        if(i == gA) stopA <- gA
      } else stopA <- 0
      text(tick.adj[c((gA>0):stopA,stopB:((gA+gB)*(gB>0)))]+0.0005*m,.26,
           gene.names.sorted[c.vec.sorted==1][c((gA>0):stopA,stopB:((gA+gB)*(gB>0)))],
           srt=90,cex=.5,font=1+2*italic)
    }
    par(adj=0.5,plt=c(.1,.9,.5,1))
    plot(c(1,m),c(0,1),xlim=c(1,m),xaxs="i",yaxs="i",xlab=x.label,ylab="",main="",col=0)
    blocks <- c(sum(local.stats<limits[1]),sum(local.stats<limits[2])+1)
    if(limits[1]> -Inf) polygon(c(m*0.0025,blocks[1],blocks[1],m*0.0025), 
                                c(0.004,0.004,0.996,0.996), border=FALSE,col=8)
    if(limits[2]< Inf) polygon(c(blocks[2],m*0.9975,m*0.9975,blocks[2]), 
                               c(0.004,0.004,0.996,0.996), border=FALSE,col=8)
    lines(c(1,m),c(0,1),lty=2, lwd=0.5)
    lines(1:m,cdf,type="s", lty=1,lwd=0.5)  
    if (!is.null(p.val) & !is.null(cat.name)){ 
      if (mean(cdf)<=0.5){
        text(m/25,0.87,cat.name,cex=0.8,font=2,adj=0)
        text(m/25,0.80,paste("p =",round(p.val,4)),cex=0.8,font=2,adj=0)
        if(substr(cat.name,1,3)=="GO:"){
          text(m/25,0.75,mget(cat.name, GOTERM)[[1]]@Term,cex=0.6,font=2,adj=0)
          text(m/25,0.70,paste("Ont:",mget(cat.name, GOTERM)[[1]]@Ontology),
               cex=0.6,font=2,adj=0)
       }
      } else {
        if(substr(cat.name,1,3)=="GO:"){
          text(24*m/25,0.27,cat.name,cex=0.8,font=2,adj=1)
          text(24*m/25,0.2,paste("p =",round(p.val,4)),cex=0.8,font=2,adj=1)
          text(24*m/25,0.15,mget(cat.name, GOTERM)[[1]]@Term,cex=0.6,font=2,adj=1)
          text(24*m/25,0.10,paste("Ont:",mget(cat.name, GOTERM)[[1]]@Ontology),
               cex=0.6,font=2,adj=1)
       } else {
          text(24*m/25,0.2,cat.name,cex=0.8,font=2,adj=1)
          text(24*m/25,0.13,paste("p =",round(p.val,4)),cex=0.8,font=2,adj=1)
       }
   }
    }
}

