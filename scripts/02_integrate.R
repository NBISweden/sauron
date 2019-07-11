#!/usr/bin/env Rscript


#############################
### LOAD/INSTALL OPTPARSE ###
#############################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')};
library(optparse)
#---------



##################################
### DEFINE PATH TO LOCAL FILES ###
##################################
cat("\nRunnign DIMENSIONLITY REDUCTION AND CLUSTERING with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object"),
  make_option(c("-c", "--columns_metadata"),      type = "character",   metavar="character",   default='none',  help="Column names in the Metadata matrix (only factors allowed, not continuous variables)"),
  make_option(c("-r", "--regress"),               type = "character",   metavar="character",   default='none',  help="Variables to be regressed out using linear modeling."),
  make_option(c("-b", "--integration_method"),    type = "character",   metavar="character",   default='none',  help="Integration method to be used. 'CCA', MNN', 'Scale' and 'Combat' are available at the moment. The batches (column names in the metadata matrix) to be removed should be provided as arguments comma separated. E.g.: 'Combat,sampling_day'. For MNN, an additional integer parameter is supplied as the k-nearest neighbour."),
  make_option(c("-p", "--PCs_use"),               type = "character",   metavar="character",   default='top,5', help="Method and threshold level for selection of significant principal components. The method should be separated from the threshold via a comma. 'top,5' will use the top 5 PCs, which is the default. 'var,1' will use all PCs with variance above 1%."),
  make_option(c("-v", "--var_genes"),             type = "character",   metavar="character",   default='Seurat,1.5',  help="Whether use 'Seurat' or the 'Scran' method for variable genes identification. An additional value can be placed after a comma to define the level of dispersion wanted for variable gene selection. 'Seurat,2' will use the threshold 2 for gene dispersions. Defult is 'Seurat,1.5'. For Scran, the user should inpup the level of biological variance 'Scran,0.2'. An additional blocking parameter (a column from the metadata) can ba supplied to 'Scran' method block variation comming from uninteresting factors, which can be parsed as 'Scran,0.2,Batch'."),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output directory")
) 
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)

col_scale <- c("grey85","navy")
if(!dir.exists(paste0(opt$output_path,"/tSNE_plots"))){dir.create(paste0(opt$output_path,"/tSNE_plots"),recursive = T)}
if(!dir.exists(paste0(opt$output_path,"/PCA_plots"))){dir.create(paste0(opt$output_path,"/PCA_plots"),recursive = T)}
#---------



##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")
initial.options <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(sub("--file=","",initial.options[grep("--file=",initial.options)]))
source( paste0(script_path,"/inst_packages.R") )
pkgs <- c("Seurat","rafalib","scran","biomaRt","scater","dplyr","RColorBrewer","dbscan","flowPeaks","scales","igraph","sva")
inst_packages(pkgs)
#---------



#############################
### LOAD Seurat.v3 OBJECT ###
#############################
cat("\nLoading/ data and metadata ...\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------



###################################
### SELECT CELLS FROM A CLUSTER ###
###################################
if (length(unlist(strsplit(opt$cluster_use,","))) >= 2 ){
  clustering_use <- as.character(unlist(strsplit(opt$cluster_use,",")))[1]
  clusters_to_select <- as.character(unlist(strsplit(opt$cluster_use,",")))[-1]
  
  if(!(clustering_use %in% colnames(DATA@meta.data))){
    cat("\nThe Clustering level was not found in the Metadata...\n")
    cat("\nAll cells will be used ...\n")
    cells_use <- rownames(DATA@meta.data)
    
  } else if(sum(clusters_to_select %in% unique(DATA@meta.data[,clustering_use])) == 0){
    cat("\nThe Cluster specifed was not found ...\n")
    cat("\nAll cells will be used ...\n")
    cells_use <- rownames(DATA@meta.data)
    
  } else {
    cells_use <- rownames(DATA@meta.data)[factor(DATA@meta.data[,clustering_use]) %in% clusters_to_select]   #Filter out cells with no assigned clusters
    
  } } else {
    cat("\nThe name of the cluster or the cluster name were not found in your data. All cells will be used ...\n")
    cells_use <- rownames(DATA@meta.data)
  }
sel <- rowSums(as.matrix(DATA@raw.data) >= 1) >= 1
DATA <- CreateSeuratObject(as.matrix(DATA@raw.data[sel,cells_use]), meta.data = DATA@meta.data[cells_use,])
#---------



#######################################
### INTEGRATE DATASETS USING COMBAT ###
#######################################
integration_method <- unlist(strsplit(opt$integration_method,","))
print(integration_method)

if ((length(integration_method) >= 2) & (casefold(integration_method[1]) == "combat") ){
  cat("\nRemoving bacthes from raw counts using COMBAT ...\n")
  
  #Defining batch variables
  batch <- factor(DATA@meta.data[,integration_method[2]])
  mod0 <- model.matrix(~1, data=as.data.frame(DATA@meta.data))
  
  #Transforming counts to log
  logdata <- log2(as.matrix(DATA@raw.data)[,rownames(DATA@meta.data)]+1)
  sum(rowSums(logdata) == 0)
  logdata <- logdata[rowSums(logdata) != 0,]
  
  #Applying Combat
  combat_data <- ComBat(dat=logdata, batch=batch, mod=mod0)
  
  #Transforming counts back to original scale (and removing negative values to 0)
  combat_data <- round(2^(combat_data)-1,0)
  sum(combat_data < 0)
  combat_data[combat_data < 0] <- 0
  DATA <- SetAssayData( object = pbmc_small, slot = "counts", new.data = combat_data, assay = "integrated" )
  DefaultAssay(DATA) <- "integrated"
  rm(combat_data,logdata,mod0);  invisible(gc())
}
#---------



####################################
### INTEGRATE DATASETS USING MNN ###
####################################
if ((length(integration_method) >= 2) & (casefold(integration_method[1]) == "mnn") ){
  cat("\nRemoving bacthes from raw counts using MNN ...\n")
  
  #Defining batch variables
  batch <- as.character(factor(DATA@meta.data[,integration_method[2]]))
  
  #Separating batch matricies 
  myinput <- list()
  for(i in unique(batch)){
    myinput[[i]] <- DATA@raw.data[,batch == i]
  }
  print(names(myinput))
  head(myinput)
  
  if(  is.na(integration_method[3]) ) { myinput[["k"]] <- 20
  } else { myinput[["k"]] <- as.numeric(integration_method[3]) }
  
  #Applying MNN correction on raw counts
  out <- do.call(mnnCorrect,args = myinput)
  mnn_cor <- do.call(cbind,out$corrected)
  
  #Converting MNN estimates back to raw counts using linear regression
  mods <- lapply(1:ncol(mnn_cor),function(x) lm(DATA@raw.data[,x] ~ mnn_cor[,x])$coefficients )
  coef2 <- setNames( t(as.data.frame(lll))[,2],colnames(DATA@raw.data))
  coef1 <- setNames( t(as.data.frame(lll))[,1],colnames(DATA@raw.data))
  
  DATA <- SetAssayData( object = pbmc_small, slot = "counts", new.data = t(round(t(mnn_cor) * coef2 + coef1,0)), assay = "integrated" )
  DefaultAssay(DATA) <- "integrated"
  rm(out, myinput);  invisible(gc())
}
#---------



####################################
### INTEGRATE DATASETS USING CCA ###
####################################
if ((length(integration_method) >= 1) & (casefold(integration_method[1]) == "cca") ){
cat("\nIntegrating datasets with CCA ...\n")
if( as.logical(opt$integrate) & (length(datasets) > 1) ){
  DATA.list <- SplitObject(DATA, split.by = "orig.ident")
  for (i in 1:length(DATA.list)) {
    DATA.list[[i]] <- NormalizeData(DATA.list[[i]], verbose = FALSE)
    DATA.list[[i]] <- FindVariableFeatures(DATA.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    gc()
  }
  DATA.anchors <- FindIntegrationAnchors(object.list = DATA.list, dims = 1:30)
  DATA <- IntegrateData(anchorset = DATA.anchors, dims = 1:30)
  DefaultAssay(DATA) <- "integrated"
  rm(DATA.list); gc()
}
}
#---------



###########################
### FIND VARIABLE GENES ###
###########################
cat("\nNormalizing and identifying highly variable genes ...\n")
DATA <- NormalizeData(DATA)

VAR_choice <- as.character(unlist(strsplit(opt$var_genes,",")))
if(casefold(VAR_choice[1]) == "no"){
  #Skip running variable gene selection and use all
  DATA@var.genes <- rownames(DATA@assays[[DefaultAssay(DATA)]]@data)
  
} else {
  
  if( (length(VAR_choice) >=2 )  &  (casefold(VAR_choice[1]) == "seurat") ){  y_cut <- as.numeric(VAR_choice[2])
  } else {  y_cut <- 2 }
  
  perc <- rowSums(as.matrix(DATA@assays[[DefaultAssay(DATA)]]@data > 0)) / ncol(DATA@assays[[DefaultAssay(DATA)]]@data)
  perc <- names(perc[ (perc < 0.80) & (perc > 10/ncol(DATA@assays[[DefaultAssay(DATA)]]@data) ) ])
  
  
  #########################################################
  ### Running SEURAT method for variable gene selection ###
  #########################################################
  cat("\nCalculating highly variable genes with Seurat ...\n")
  #Defining the variable genes based on the mean gene expression abothe the 5% quantile and the dispersion above 2.
  DATA <- FindVariableGenes(object = DATA, mean.function = ExpMean, dispersion.function = LogVMR, y.cutoff = y_cut,num.bin = 200)
  m <- max(quantile(DATA@hvg.info$gene.mean,probs = c(.025)) , 0.01)
  DATA <- FindVariableGenes(object = DATA, mean.function = ExpMean, dispersion.function = LogVMR, y.cutoff = y_cut,num.bin = 200,x.low.cutoff = m)
  
  DATA@var.genes <- DATA@var.genes[DATA@var.genes %in% perc]
  DATA@hvg.info$use <- rownames(DATA@hvg.info) %in% DATA@var.genes
  write.csv2(DATA@hvg.info, paste0(opt$output_path,"/HVG_info_seurat.csv"))
  
  png(filename = paste0(opt$output_path,"/Var_gene_selection_seurat.png"),width = 700,height = 750,res = 150)
  plot(log2(DATA@hvg.info$gene.mean),DATA@hvg.info$gene.dispersion.scaled,cex=.1,main="HVG selection",
       col=ifelse(rownames(DATA@hvg.info)%in% DATA@var.genes,"red","black" ),ylab="scaled.dispersion",xlab="log2(avg. expression)")
  abline(v=log2(m),h=y_cut,lty=2,col="grey20",lwd=1)
  invisible(dev.off())
  #---------
  
  if( (length(VAR_choice) >=2 )  &  (casefold(VAR_choice[1]) == "scran") ){  y_cut <- as.numeric(VAR_choice[2])
  } else {  y_cut <- 0.15 }
  
  
  
  ########################################################
  ### Running SCRAN method for variable gene selection ###
  ########################################################
  cat("\nCalculating highly variable genes with Scran ...\n")
  if( (length(VAR_choice)==3) & (VAR_choice[3] %in% colnames(DATA@meta.data)) ){
    cat("\nBlocking factor detected ...\n")
    blk <- DATA@meta.data[,VAR_choice[3]]
    fit <- trendVar(DATA@data,loess.args=list(span=0.05), block=blk)
    fit$vars <- apply(fit$vars,1, function(x) {prod(x)^(1/length(x))} )
    fit$means <- apply(fit$means,1, function(x) {prod(x)^(1/length(x))} )
  } else { fit <- trendVar(DATA@data,loess.args=list(span=0.05)) }
  
  hvgs <- decomposeVar(DATA@data, fit)
  hvgs <- as.data.frame(hvgs[order(hvgs$bio, decreasing=TRUE),])
  
  #The minimum variance. The minimum variance needs to be so that at least 1% of the cells express that gene
  n <- ncol(DATA@data)
  min_var <- var(sample(size = n,x = c(1,0),replace = T,prob = c((n/100)/n,(n - (n/100))/n) ))
  myvars <- rownames(hvgs)[ (fit$vars > min_var) & (hvgs$bio > y_cut) & (hvgs$FDR < 0.01 ) ]
  myvars <- myvars[myvars %in% perc]
  hvgs$use <- rownames(hvgs) %in% myvars
  write.csv2(hvgs, paste0(opt$output_path,"/HVG_info_scran.csv"))
  
  png(filename = paste0(opt$output_path,"/Var_fit_scran.png"),width = 700,height = 750,res = 150)
  TF <- names(fit$means) %in% myvars
  plot( c(fit$means), c(fit$vars) ,xlab="mean",ylab="biological variance",pch=16,
        col=ifelse(TF ,"red","grey30"),
        cex=ifelse(TF ,.5,.2) , main="SCRAN")
  curve(fit$trend(x), col="red", lwd=2, add=TRUE)
  invisible(dev.off())
  
  png(filename = paste0(opt$output_path,"/Var_genes_scran.png"),width = 700,height = 750,res = 150)
  plot( log2(hvgs$mean) , log2(hvgs$bio+1) ,xlab="log2(mean)",ylab="log2(bio.var+1)",pch=16,
        col=ifelse(hvgs$use,"red","grey30"),ylim=c(-0.1,2),
        cex=ifelse(hvgs$use,.5,.2) ,main="SCRAN")
  invisible(dev.off())
  #---------
  
  
  if(casefold(VAR_choice[1]) == "scran"){
    cat("\nSCRAN was the method chosen for downstream procedure ...\n")
    DATA@hvg.info <- hvgs
    DATA@var.genes <- myvars
  } else {
    cat("\nSEURAT was the method chosen for downstream procedure ...\n")
  }
}
#---------





#############################################
### Scaling data and regressing variables ###
#############################################
cat("\nScaling data and regressing uninteresting factors ...\n")
integration_method <- unlist(strsplit(opt$integration_method,","))
vars <- as.character(unlist(strsplit(opt$regress,",")))

if ((length(integration_method) >= 2) & (integration_method[1] == "Scale") ){
  vars <- unique(c(vars, integration_method[2:length(integration_method)]))
}

DATA <- ScaleData(DATA,vars.to.regress = vars)
#---------



###################################
### SAVING RAW Seurat.v3 OBJECT ###
###################################
cat("\nSaving the RAW Seurat object ...\n")
saveRDS(DATA, file = paste0(opt$output_path,"/Seurat_object.rds") )
#---------



#############################
### SYSTEM & SESSION INFO ###
#############################
#---------
print_session_info()
#---------

