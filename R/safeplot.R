safeplot <-
function(safe = NULL, cat.name = "", limits = NULL,
                       c.vec = NULL, local.stats = NULL, gene.names = NULL,
                       rank = TRUE, x.limits = NULL, c.thresh = 0,
                       colors = NULL, x.ticks = NULL, t.cex = 1,
                       p.val = NULL, cat.desc = NULL, title = "",...){
   if(!is.null(safe)){
      local.stats <- safe@local.stat
      C.names <- names(safe@global.stat)
      if(cat.name=="") cat.name <- C.names[order(safe@global.pval)][1]
      if(sum(is.null(limits))){
         limits <- c(0,Inf)
         if(min(local.stats)<0)  limits[1] <- min(local.stats[(safe@local.pval > 0.05)&(local.stats < 0)])
         limits[2] <- max(local.stats[(safe@local.pval > 0.05)&(local.stats > 0)])
         cat(paste("Shaded limits are the largest local statistics with p > 0.05\n",
                   "  Limits = (",round(limits[1],3),",",round(limits[2],3),")\n"))
      }
      c.vec <- as.matrix(safe@C.mat[,C.names == cat.name])
      if(is.null(p.val)) p.val <- sigfig(safe@global.pval[C.names == cat.name],4)
   }
   if(is.null(p.val)) p.val <- NA
   if(is.null(colors)) colors <- rep(1,length(local.stats))
   if(is.null(limits)) limits <- range(local.stats)
   if(is.null(gene.names)) gene.names <- names(local.stats)
   if(is.null(cat.desc)){
       if(substr(cat.name,1,5) == "KEGG:") {
           require(KEGG.db); cat.desc <-
             mget(substr(cat.name,6,10),KEGGPATHID2NAME)[[1]]}
       if(substr(cat.name,1,3) == "GO:") {
           require(GO.db); cat.desc <-
             Term(mget(cat.name,GOTERM)[[1]])}
       if(substr(cat.name,1,2) == "PF"){
           require(PFAM.db); cat.desc <-
             mget(cat.name,PFAMID,ifnotfound=NA)[[1]]}
   }
   g <- sum(c.vec);m <- length(c.vec)
   if(rank){
       x <- 1:m
       ecdf <-  approx(x,(1:m)/m)
   } else {
       if(is.null(x.limits)) x.limits <- range(local.stats)
       ecdf <- approx(sort(local.stats),(1:m)/m)
       x <- sapply(local.stats,min,x.limits[2])
       x <- sort(sapply(x,max,x.limits[1]))
   }

   c.vec.sorted <- c.vec[order(local.stats)]
   gene.names.sorted <- gene.names[order(local.stats)]
   colors.sorted <- colors[order(local.stats)]
   cdf <- cumsum(c.vec.sorted)/ g

   if(prod(c.vec.sorted %in% 0:1)) c.vec.sorted <- c.vec.sorted * 0.25

   nf <- layout(matrix(c(1:3),3,1,byrow=TRUE), c(4,4,4), c(0.5,1,2), TRUE)
   par(plt=c(.1,.9,0,1),adj=0)
   plot(0,0,axes=FALSE,ylab="",xlab= "",col=0);text(0,0,title,adj=0.5,cex=2)
   plot(x,rep(0.5,m),ylim=c(0,1),xaxs="i",axes=FALSE,ylab="",xlab= "",col=0,yaxs="i")
   for(i in 1:m) if(c.vec.sorted[i]>0) points(rep(x[i],2),c(0,0.2*c.vec.sorted[i]),
                                       type="l",lwd=1.5)

      if(rank) zero   <- sum(local.stats<0) + 0.5 else zero <- 0
      tick   <- x[c.vec.sorted>c.thresh]
      tick.c <- colors.sorted[c.vec.sorted>c.thresh]
      tick.x <- c.vec.sorted[c.vec.sorted>c.thresh]
      gA <- sum(tick<zero)
      gB <- sum(tick>zero)
      step <- 0.02 * diff(range(x))
      tick.adj <- c(rep(max(tick[1],min(x)+step/2),gA),rep(min(tick[gA+gB],max(x)-step/2),gB))
      i <- 1;  stopA <- 0
      if(gA>0){
        tick.adj[i] <- min(tick.adj[i],zero-step/2)
        lines(c(tick.adj[i],tick[i]), c(0.25, 0.2*tick.x[i]),col=tick.c[i])
        if(gA==1) stopA <- 1 else while(stopA == 0){
          i <- i + 1
          tick.adj[i]<- max(tick.adj[i-1]+step,tick[i])
          if(tick.adj[i]<(zero-step/2)){
            lines(c(tick.adj[i],tick[i]), c(0.25, 0.2*tick.x[i]),lwd=0.25,col=tick.c[i])
            if(i == gA) stopA <- gA
          } else stopA <- i-1
        }
      }
      i <- gA + gB;  stopB <- 0
      if(gB>0){
         tick.adj[i] <- max(tick.adj[i],zero+step/2)
         lines(c(tick.adj[i],tick[i]), c(0.25, 0.2*tick.x[i]),col=tick.c[i])
         if(gB==1) stopB <- gA+1 else while(stopB == 0){
           i <- i - 1
           tick.adj[i]<-min(tick.adj[i+1]-step,tick[i])
           if(tick.adj[i]>(zero+step/2)){
             lines(c(tick.adj[i],tick[i]), c(0.25,0.2*tick.x[i] ),lwd=0.25,col=tick.c[i])
             if(i == gA+1) stopB <- gA+1
           } else stopB <- i+1
         }
      }
      text(tick.adj[c((gA>0):stopA,stopB:((gA+gB)*(gB>0)))]+0.0005*max(x),.26,
           gene.names.sorted[c.vec.sorted>c.thresh][c((gA>0):stopA,stopB:((gA+gB)*(gB>0)))],
           col= colors.sorted[c.vec.sorted>c.thresh][c((gA>0):stopA,stopB:((gA+gB)*(gB>0)))],
           srt=90,cex=t.cex,font=1)
    par(adj=0.5,plt=c(.1,.9,.5,1),cex.lab=1.5)
    if(rank) xl <- "Local statistics (on ranked scale)" else xl <- "Local statistics"
    plot(range(x),c(0,1),xlim=range(x),xaxs="i",yaxs="i",axes=F,xlab=xl,ylab="",main="",col=0)

    if(rank){
      blocks <- c(sum(local.stats<limits[1]),sum(local.stats<limits[2])+1)
    } else blocks <- limits
    if(limits[1]> -Inf) polygon(c(x[1],blocks[1],blocks[1],x[1]),
                                c(0.004,0.004,0.996,0.996), border=FALSE,col=8)
    if(limits[2]< Inf) polygon(c(blocks[2],x[m],x[m],blocks[2]),
                               c(0.004,0.004,0.996,0.996), border=FALSE,col=8)
    axis(1);axis(2,seq(0,1,0.25),c("0","0.25","0.5","0.75","1"));box()

    lines(ecdf,lty=2, lwd=0.5)
    if(rank & !is.null(x.ticks)){
       x.n <- sapply(x.ticks,function(x) sum(local.stats < x))
       for(i in 1:length(x.n)) {
        lines(rep(x.n[i],2),0:1,lwd=0.5,lty=2)
        text(x.n[i],0,paste("t=",x.ticks[i],sep=""),adj=c(-0.2,-0.2),cex=0.75)
       }
   }
    lines(ecdf,lty=2, lwd=0.5)

    d <- duplicated(cdf)
    lines(c(min(x),x[!d],max(x)),c(0,cdf[!d],0.98),type="s", lty=1,lwd=2)
    if (!is.null(p.val) & !is.null(cat.name)){
      if (mean(cdf)<=0.5){
        text(min(x)+step,0.90,cat.name,cex=1.5,font=2,adj=c(0,0.5))
        if(!is.na(p.val)) text(min(x)+step,0.80,paste("p =",p.val),cex=1.25,font=2,adj=c(0,0.5))
        if(!is.null(cat.desc)) text(min(x)+step,0.71,cat.desc,cex=1,font=2,adj=c(0,0.5))
     } else {
        text(max(x)-step,0.27,cat.name,cex=1.5,font=2,adj=c(1,0.5))
        if(!is.na(p.val)) text(max(x)-step,0.17,paste("p =",p.val),cex=1.25,font=2,adj=c(1,0.5))
        if(!is.null(cat.desc)) text(max(x)-step,0.08,cat.desc,cex=1,font=2,adj=c(1,0.5))
    }
  }
}
