#!/usr/bin/env Rscript

compute_hvgs <- function(DATA,VAR_choice,output_path){
  if(!dir.exists(output_path)){dir.create(output_path,recursive = T)}
  
  VAR_choice <- as.character(unlist(strsplit(opt$var_genes,",")))
  if(casefold(VAR_choice[1]) == "no"){
    #Skip running variable gene selection and use all
    DATA@assays[[DefaultAssay(DATA)]]@var.features <- rownames(DATA@assays[[DefaultAssay(DATA)]]@data)
    
  } else {
    
    perc1 <- rowSums(as.matrix(DATA@assays[[DefaultAssay(DATA)]]@data > 0)) / ncol(DATA@assays[[DefaultAssay(DATA)]]@data)
    perc <- names(perc1)[ (perc1 < 0.90) & (perc1 > 10/ncol(DATA@assays[[DefaultAssay(DATA)]]@data) ) ]
    
    #########################################################
    ### Running SEURAT method for variable gene selection ###
    #########################################################
    cat("\nCalculating highly variable genes with Seurat ...\n")
    if( (length(VAR_choice) >=2 )  &  (casefold(VAR_choice[1]) == "seurat") ){  y_cut <- as.numeric(VAR_choice[2])
    } else {  y_cut <- 1 }
    
    #Defining the variable genes based on the mean gene expression abothe the 5% quantile and the dispersion above 2.
    DATA <- FindVariableFeatures(object = DATA)
    m <- max(quantile(DATA@assays[[DefaultAssay(DATA)]]@meta.features$vst.mean,probs = c(.025)) , 0.01)
    
    DATA@assays[[DefaultAssay(DATA)]]@var.features <- DATA@assays[[DefaultAssay(DATA)]]@var.features[DATA@assays[[DefaultAssay(DATA)]]@var.features %in% perc]
    DATA@assays[[DefaultAssay(DATA)]]@meta.features$use <- rownames(DATA@assays[[DefaultAssay(DATA)]]@meta.features) %in% DATA@assays[[DefaultAssay(DATA)]]@var.features
    write.csv2(DATA@assays[[DefaultAssay(DATA)]]@meta.features, paste0(output_path,"/HVG_info_seurat.v3.csv"))
    
    png(filename = paste0(output_path,"/Var_vst_exp_disp_gene_selection_seurat.v3.png"),width = 700,height = 750,res = 150)
    plot(log2(DATA@assays[[DefaultAssay(DATA)]]@meta.features$vst.variance.expected),DATA@assays[[DefaultAssay(DATA)]]@meta.features$vst.variance.standardized,cex=.1,main="HVG selection",
         col=ifelse(rownames(DATA@assays[[DefaultAssay(DATA)]]@meta.features)%in% DATA@assays[[DefaultAssay(DATA)]]@var.features,"red","black" ),ylab="vst.variable",xlab="log2(var.expected)")
    abline(v=log2(m),h=y_cut,lty=2,col="grey20",lwd=1)
    invisible(dev.off())
    
    png(filename = paste0(output_path,"/Var_vst_mean_disp_gene_selection_seurat.v3.png"),width = 700,height = 750,res = 150)
    plot(log2(DATA@assays[[DefaultAssay(DATA)]]@meta.features$vst.mean),DATA@assays[[DefaultAssay(DATA)]]@meta.features$vst.variance.standardized,cex=.1,main="HVG selection",
         col=ifelse(rownames(DATA@assays[[DefaultAssay(DATA)]]@meta.features)%in% DATA@assays[[DefaultAssay(DATA)]]@var.features,"red","black" ),ylab="vst.variable",xlab="log2(mean expression)")
    abline(v=log2(m),h=y_cut,lty=2,col="grey20",lwd=1)
    invisible(dev.off())
    #---------
    
    
    
    ########################################################
    ### Running SCRAN method for variable gene selection ###
    ########################################################
    cat("\nCalculating highly variable genes with Scran ...\n")
    
    
    if( (length(VAR_choice)==3) & (VAR_choice[3] %in% colnames(DATA@meta.data)) ){
      cat("\nBlocking factor detected ...\n")
      blk <- DATA@meta.data[,VAR_choice[3]]
      fit <- trendVar(DATA@assays[[DefaultAssay(DATA)]]@data,loess.args=list(span=0.05), block=blk)
      fit$vars <- apply(fit$vars,1, function(x) {prod(x)^(1/length(x))} )
      fit$means <- apply(fit$means,1, function(x) {prod(x)^(1/length(x))} )
    } else { fit <- trendVar(DATA@assays[[DefaultAssay(DATA)]]@data,loess.args=list(span=0.05)) }
    
    hvgs <- decomposeVar(DATA@assays[[DefaultAssay(DATA)]]@data, fit)
    hvgs <- as.data.frame(hvgs[order(hvgs$bio, decreasing=TRUE),])
    hvgs$percentage <- perc1
    
    if( (length(VAR_choice) >=2 )  &  (casefold(VAR_choice[1]) == "scran") ){  y_cut <- as.numeric(VAR_choice[2])
    } else {  y_cut <- quantile(log2(hvgs$bio+1),probs=.9,na.rm=T) }
    ?quantile
    
    #The minimum variance. The minimum variance needs to be so that at least 10 of the cells express that gene
    n <- ncol(DATA@assays[[DefaultAssay(DATA)]]@data)
    min_var <- var(c(rep(1, 10 ),rep(0,round( n - (n/100)) )))
    myvars <- rownames(hvgs)[ (fit$vars > min_var) & (hvgs$bio > y_cut) & (hvgs$FDR < 0.01 ) ]
    myvars <- myvars[myvars %in% perc]
    hvgs$use <- rownames(hvgs) %in% myvars
    write.csv2(hvgs, paste0(output_path,"/HVG_info_scran.csv"))
    
    png(filename = paste0(output_path,"/Var_fit_scran.png"),width = 2*700,height = 750,res = 150)
    par(mfrow=c(1,2))
    TF <- names(fit$means) %in% myvars
    plot( c(fit$means), c(fit$vars) ,xlab="mean",ylab="biological variance",pch=16,
          col=ifelse(TF ,"red","grey30"),
          cex=ifelse(TF ,.5,.2) , main="SCRAN")
    curve(fit$trend(x), col="red", lwd=2, add=TRUE)
    plot( log2(c(fit$means)), c(fit$vars) ,xlab="mean",ylab="biological variance",pch=16,
          col=ifelse(TF ,"red","grey30"),
          cex=ifelse(TF ,.5,.2) , main="SCRAN")
    invisible(dev.off())
    
    png(filename = paste0(output_path,"/Var_genes_scran_percent.png"),width = 3*700,height = 750,res = 150)
    par(mfrow=c(1,3))
    plot( log2(hvgs$mean) , log2(hvgs$bio+1) ,xlab="log2(mean)",ylab="log2(bio.var+1)",pch=16,
          col=ifelse(hvgs$use,"red","grey30"),ylim=c(-0.1,2),
          cex=ifelse(hvgs$use,.5,.2) ,main="SCRAN")
    plot( hvgs$percentage , log2(hvgs$bio+1) ,xlab="detection rate",ylab="log2(bio.var+1)",pch=16,
          col=ifelse(hvgs$use,"red","grey30"),ylim=c(-0.1,2),
          cex=ifelse(hvgs$use,.5,.2) ,main="SCRAN")
    plot( log2(hvgs$mean), hvgs$percentage,xlab="log2(mean)",ylab="detection rate",pch=16,
          col=ifelse(hvgs$use,"red","grey30"),ylim=c(-0.1,1),
          cex=ifelse(hvgs$use,.5,.2) ,main="SCRAN")
    invisible(dev.off())
    #---------
    
    
    if(casefold(VAR_choice[1]) == "scran"){
      cat("\nSCRAN was the method chosen for downstream procedure ...\n")
      DATA@assays[[DefaultAssay(DATA)]]@meta.features <- hvgs
      DATA@assays[[DefaultAssay(DATA)]]@var.features <- myvars
    } else {
      cat("\nSEURAT was the method chosen for downstream procedure ...\n")
    }
  }
  return(DATA)
}