library(scales)
#MAIN Violin plot function
#---------------
violins <- function(data, gene, clustering, plot_points=T,plot_y_axis=T,plot_x_axis=T,smooth=2,method="log",points_method="log",col="default",
                    pt.col="grey",bw=.7,max_points=200,assay="RNA",ylab="expression",cex.main=1.7,main=gene,...){
  #par(mar=c(2,3,2,1))
  n <- length(unique(data@meta.data[,clustering]))
  my_max <- max(data@assays[[assay]]@data[gene,],1)*1.1
  plot(c(.4,n+.6),c(-1,-1), ylim=c(-.1,my_max),...,ylab="",type="n" ,frame.plot = F,yaxs="i",xaxs="i",las=1,xlab="",main=main,xaxt = "n",yaxt = "n",cex.main=cex.main)
  mtext(side = 2, text = ylab, line = 2,las=3)
  #col_pal <- hue_pal()(length(unique(data@meta.data[,clustering])))
  col_pal <- c(scales::hue_pal()(8),RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(8,"Set2") )
  col_pal <- col_pal[sort(unique(data@meta.data[,clustering]))]
  
  if(col=="default"){col <- paste0(col_pal,95)} else { col <- rep(col, length(unique(data@meta.data[,clustering])) )}
  
  for(i in 1:length(unique(data@meta.data[,clustering]))){
    cl <- sort(unique(data@meta.data[,clustering]))[i]
    #if(plot_points){
      #points(runif(sum(DATA@meta.data[,clustering] == cl),min = i-.4,max = i+.4),DATA@data[gene,DATA@meta.data[,clustering] == cl],cex=.3,pch=16,col="grey40")
      #points(rnorm(sum(DATA@meta.data[,clustering] == cl),mean = i,sd = .12),DATA@data[gene,DATA@meta.data[,clustering] == cl],cex=.3,pch=16,col="grey60")
    #}
    x <- data@assays[[assay]]@data[gene,data@meta.data[,clustering] == cl]
    suppressWarnings(suppressMessages( try(draw_violin(x,at = i,col = col[i], smooth=smooth,plot_points=plot_points,method=method,points_method="proportional",
               bw = bw,border =  "grey20",max_points=max_points)) ))
    #paste0(col_pal[i])
    #vioplot(x,at = i,add=T,col = paste0(col_pal[i],95),
            #drawRect = F,wex = 1,h = .01, border =  paste0(col_pal[i]))
  }
  lines(c(0,n+.6),c(my_max,my_max),col="white",lwd=6,xpd=T)
  lines(c(n+.6,n+.6),c(0,my_max),col="white",lwd=6,xpd=T)
  abline(h=-.1,v=.4,xpd=F,lwd=2)
  if(plot_x_axis){
    axis(1, at=1:n, labels=sort(unique(data@meta.data[,clustering])),cex.axis=1.5)
  }
  if(plot_y_axis){
      axis(2, at=1:10, labels=1:10,cex.axis=1.5,las=1)
  }
}
#---------------


#Function to calculate violin density
#---------------
draw_violin <- function(x,method="log",plot_points=F,points_method="proportional",smooth=2,col="grey",border="grey",at=1,pt.col="grey",bw=0.45,max_points=200){
  r <- sum(x!=0)/length(x)
  if(plot_points){
    if(points_method == "proportional"){points(rnorm(length(x),mean = at,r/5),x,cex=.5,col=pt.col,pch=16)}
    if(points_method == "uniform"){points(rnorm(length(x),mean = at,sd = .12),x,cex=.5,col=pt.col,pch=16)}
  }

  if(method == "uniform"){
    a <- density(x,bw = .25*smooth,n = 200)
    ys <- a$x
    xs <- a$y/max(a$y)*bw/2
  }
  
  if(method == "log"){
    x2 <- log2(log2(x+1)+1)
    a <- density(x2,bw = .05*smooth,n = 200)
    ys <- 2^(2^a$x-1)-1
    xs <- a$y/max(a$y)*bw/2
  }
  
  if(method == "mixed"){
    d1 <- density(x[x==0],bw = .01*smooth,n = 200)
    d1$y <- d1$y/max(d1$y)*(1-r)
    d2 <- density(c(x[x!=0]),bw = .2,n = 200)
    d2$y <- d2$y/max(d2$y)*(r)
    ys <- c(d1$x,d2$x)
    xs <- c(d1$y,d2$y)*bw
  }
  ulim <- max(quantile(x,.99),0.1)
  llim <- max(min(x),-0.1)
  polygon( c(xs[ys<ulim & ys>llim] , -rev(xs[ys<ulim& ys>llim]) )+at, 
           c(ys[ys<ulim & ys>llim], rev(ys[ys<ulim& ys>llim])) ,col = col,border = border)
}
#---------------
