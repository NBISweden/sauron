#!/usr/bin/env Rscript

### LOAD OPTPARSE
#---------
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')};
library(optparse)
#---------


### DEFINE PATH TO LOCAL FILES
#---------
cat("\nRunnign DIMENS> REDUCTION AND CLUSTERING with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object"),
  make_option(c("-c", "--columns_metadata"),      type = "character",   metavar="character",   default='none',  help="Column names in the Metadata matrix (only factors allowed, not continuous variables)"),
  make_option(c("-r", "--regress"),               type = "character",   metavar="character",   default='none',  help="Variables to be regressed out using linear modeling."),
  make_option(c("-r", "--batch_method"),          type = "character",   metavar="character",   default='none',  help="Batch-correction method to be used. 'MNN', 'Scale' and 'Combat' are available at the moment. The batches (column names in the metadata matrix) to be removed should be provided as arguments comma separated. E.g.: 'Combat,sampling_day'. For MNN, an additional integer parameter is supplied as the k-nearest neighbour."),
  make_option(c("-p", "--PCs_use"),               type = "character",   metavar="character",   default='top,5', help="Method and threshold level for selection of significant principal components. The method should be separated from the threshold via a comma. 'top,5' will use the top 5 PCs, which is the default. 'var,1' will use all PCs with variance above 1%."),
  make_option(c("-v", "--var_genes"),             type = "character",   metavar="character",   default='Seurat,1.5',  help="Whether use 'Seurat' or the 'Scran' method for variable genes identification. An additional value can be placed after a comma to define the level of dispersion wanted for variable gene selection. 'Seurat,2' will use the threshold 2 for gene dispersions. Defult is 'Seurat,1.5'. For Scran, the user should inpup the level of biological variance 'Scran,0.2'. An additional blocking parameter (a column from the metadata) can ba supplied to 'Scran' method block variation comming from uninteresting factors, which can be parsed as 'Scran,0.2,Batch'."),
  make_option(c("-s", "--cluster_use"),           type = "character",   metavar="character",   default='none',  help="The clustering method and cluster to select for analysis"),
  make_option(c("-f", "--aux_functions_path"),    type = "character",   metavar="character",   default='none',  help="File with supplementary functions"),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output directory")
) 
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)
#---------



### DEFINE PATH TO LOCAL FILES
#---------
col_scale <- c("grey85","navy")
if(!dir.exists(paste0(opt$output_path,"/tSNE_plots"))){dir.create(paste0(opt$output_path,"/tSNE_plots"),recursive = T)}
if(!dir.exists(paste0(opt$output_path,"/PCA_plots"))){dir.create(paste0(opt$output_path,"/PCA_plots"),recursive = T)}
#---------



### LOAD LIBRARIES
#---------
cat("\nLoading/installing libraries ...\n")
source(opt$aux_functions_path)
pkgs <- c("Seurat","rafalib","scran","biomaRt","scater","dplyr","RColorBrewer","dbscan","flowPeaks","scales","igraph","sva")
inst_packages(pkgs)
#---------




### LOAD Seurat OBJECT 
#---------
cat("\nLoading/ data and metadata ...\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------




### Filter cells from a specific cluster
#---------
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
#---------





### Removing batch effects using ComBat from SVA package (on raw counts)
#---------
batch_method <- unlist(strsplit(opt$batch_method,","))

if ((length(batch_method) >= 2) & (batch_method[1] == "Combat") ){
  cat("\nRemoving bacthes from raw counts using ComBat ...\n")
  
  #Defining batch variables
  batch <- factor(DATA@meta.data[,batch_method[2]])
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
  DATA@raw.data <- combat_data
  rm(combat_data,logdata,mod0)
}
#---------






### Remove batch effects using MNN (on raw counts)
#---------
batch_method <- unlist(strsplit(opt$batch_method,","))

if ((length(batch_method) >= 3) & (batch_method[1] == "MNN") ){
  cat("\nRemoving bacthes from raw counts using MNN ...\n")

  #Defining batch variables
  batch <- factor(DATA@meta.data[,batch_method[2]])
  
  #Separating batch matricies 
  myinput <- list()
  for(i in unique(batch)){
    myinput[[i]] <- DATA@raw.data[,batch == i]
  }
  myinput[["k"]] <- as.numeric(batch_method[3])
  
  #Applying MNN correction on raw counts
  out <- do.call(mnnCorrect,args = myinput)
  mnn_cor <- do.call(cbind,out$corrected)
  
  #Converting MNN estimates back to raw counts using linear regression
  mods <- lapply(1:ncol(mnn_cor),function(x) lm(DATA@raw.data[,x] ~ mnn_cor[,x])$coefficients )
  coef2 <- setNames( t(as.data.frame(lll))[,2],colnames(DATA@raw.data))
  coef1 <- setNames( t(as.data.frame(lll))[,1],colnames(DATA@raw.data))
  
  DATA@raw.data <- t(round(t(mnn_cor) * coef2 + coef1,0))
  rm(out, myinput)
}
#---------





### Normalizing and finding highly variable genes
#---------
cat("\nNormalizing and identifying highly variable genes ...\n")
sel <- rowSums(as.matrix(DATA@raw.data) >= 2) >= 5
DATA <- CreateSeuratObject(as.matrix(DATA@raw.data[sel,cells_use]), meta.data = DATA@meta.data[cells_use,])
DATA <- NormalizeData(DATA)

VAR_choice <- as.character(unlist(strsplit(opt$var_genes,",")))
if(VAR_choice[1] == "no"){
  #Skip running variable gene selection and use all
  DATA@var.genes <- rownames(DATA@data)
  
} else {
  if(VAR_choice[1] == "Scran"){
    #Running SCRAN method for variable gene selection
    y_cut <- as.numeric(VAR_choice[2])

    if( (length(VAR_choice)==3) ){
      blk <- DATA@meta.data[,VAR_choice[3]]
      fit <- trendVar(DATA@data,loess.args=list(span=0.05), block=blk)
    } else { fit <- trendVar(DATA@data,loess.args=list(span=0.05)) }
    
    plot(fit$mean,fit$var,cex=0.2)
    curve(fit$trend(x), col="red", lwd=2, add=TRUE)
    
    hvgs <- decomposeVar(DATA@data, fit)
    hvgs <- cbind(hvgs,perc_bio= hvgs$bio / hvgs$total)
    hvgs <- as.data.frame(hvgs[order(hvgs$bio, decreasing=TRUE),])
    
    DATA@hvg.info <- hvgs
    DATA@var.genes <- rownames(hvgs)[(hvgs$FDR < 0.001) & (hvgs$bio > y_cut)]

  } else {
    #Running SEURAT method for variable gene selection
    y_cut <- as.numeric(VAR_choice[2])
    #Defining the variable genes based on the mean gene expression abothe the 5% quantile and the dispersion above 2.
    DATA <- FindVariableGenes(object = DATA, mean.function = ExpMean, dispersion.function = LogVMR, y.cutoff = y_cut,num.bin = 200)
    m <- max(quantile(DATA@hvg.info$gene.mean,probs = c(.025)) , 0.01)
    DATA <- FindVariableGenes(object = DATA, mean.function = ExpMean, dispersion.function = LogVMR, y.cutoff = y_cut,num.bin = 200,x.low.cutoff = m)
    
    png(filename = paste0(opt$output_path,"/Var_gene_selection.png"),width = 700,height = 750,res = 150)
    plot(log2(DATA@hvg.info$gene.mean),DATA@hvg.info$gene.dispersion.scaled,cex=.1,main="HVG selection",
         col=ifelse(rownames(DATA@hvg.info)%in% DATA@var.genes,"red","black" ),ylab="scaled.dispersion",xlab="log2(avg. expression)")
    abline(v=log2(m),h=y_cut,lty=2,col="grey20",lwd=1)
    dev.off()
  }
  write.csv2(DATA@hvg.info, paste0(opt$output_path,"/HVG_info.csv"))
}
#---------





### Scaling data and regressing variables
#---------
cat("\nScaling data and regressing uninteresting factors ...\n")
batch_method <- unlist(strsplit(opt$batch_method,","))
vars <- as.character(unlist(strsplit(opt$regress,",")))

if ((length(batch_method) >= 2) & (batch_method[1] == "Scale") ){
  vars <- unique(c(vars, batch_method[2:length(batch_method)]))
}

DATA <- ScaleData(DATA,vars.to.regress = vars)
#---------




### Running PCA
#---------
cat("\nRunning PCA ...\n")
DATA <- RunPCA(DATA, do.print = F, pcs.compute = 100)
var_expl <- (DATA@dr$pca@sdev^2)/sum(DATA@dr$pca@sdev^2)

PC_choice <- as.character(unlist(strsplit(opt$PCs_use,",")))
if(PC_choice[1] == "var"){
  top_PCs <- sum( var_expl > as.numeric(PC_choice[2])/100 )
} else if(PC_choice[1] == "top"){
  top_PCs <- as.numeric(PC_choice[2])
} 

png(filename = paste0(opt$output_path,"/PCA_plots/Varianc_explained_PC.png"),width = 1500,height =1200,res = 200)
plot( var_expl,yaxs="i",bg="grey",pch=21,type="l",ylab="% Variance",xlab="PCs",main="all cells")
points( var_expl,bg=c(rep("orange",top_PCs),rep("grey",100-top_PCs)),pch=21)
invisible(dev.off())

cat("\nRunning BH-tSNE ...\n")
if(file.exists(paste0(opt$output_path,"/tSNE_plots/tSNE_coordinates.csv"))){
  cat("\nPre-computed tSNE found and will be used:\n",paste0(opt$output_path,"/tSNE_plots/tSNE_coordinates.csv"),"\n")
  DATA <- SetDimReduction(DATA, reduction.type = "tsne",  slot = "cell.embeddings", new.data = as.matrix(read.csv2(paste0(opt$output_path,"/tSNE_plots/tSNE_coordinates.csv"),row.names = 1)) )
  DATA <- SetDimReduction(DATA, reduction.type = "tsne",  slot = "key", new.data = "tSNE_")
} else {
  cat("\nPre-computed tSNE NOT found. Computing tSNE ...\n")
  DATA <- RunTSNE(object = DATA, perplexity=30, max_iter=2000,theta=0,eta=2000,exaggeration_factor=12,dims.use = 1:top_PCs,verbose = T)
  write.csv2(DATA@dr$tsne@cell.embeddings, paste0(opt$output_path,"/tSNE_plots/tSNE_coordinates.csv"))
}
#---------




#Plotting PCA and tSNE plots for metadata
#---------
cat("\nPlotting PCA and tSNE plot for variables ...\n")

for(j in c("nUMI","nGene")){
  if(j %in% colnames(DATA@meta.data)){
    png(filename = paste0(opt$output_path,"/PCA_plots/PCA_",j,".png"),width = 600,height = 600,res = 150)
    print(FeaturePlot(object = DATA, features.plot = j, cols.use = col_scale,pt.size = .5,reduction.use = "pca"))
    invisible(dev.off())
  }}

for(j in as.character(unlist(strsplit(opt$columns_metadata,",")))){
  png(filename = paste0(opt$output_path,"/PCA_plots/PCA_",j,".png"),width = 700,height = 600,res = 150)
  print(PCAPlot(object = DATA,group.by=j,pt.size = .3))
  invisible(dev.off())
}

for(j in as.character(unlist(strsplit(opt$columns_metadata,",")))){
  png(filename = paste0(opt$output_path,"/tSNE_plots/tSNE_",j,".png"),width = 700,height = 600,res = 150)
  print(TSNEPlot(object = DATA,group.by=j,pt.size = .3))
  invisible(dev.off())
  }

for(j in c("nUMI","nGene","S.Score","G2M.Score","percent.Rpl","percent.Rps")){
  if(j %in% colnames(DATA@meta.data)){
  png(filename = paste0(opt$output_path,"/tSNE_plots/tSNE_",j,".png"),width = 600,height = 600,res = 150)
  print(FeaturePlot(object = DATA, features.plot = j, cols.use = col_scale,pt.size = .5))
  invisible(dev.off())
  }}

png(filename = paste0(opt$output_path,"/tSNE_plots/tSNE_density_auto.png"),width = 600,height =650,res = 100)
plot(DATA@dr$tsne@cell.embeddings,pch=16,cex=0.5, col=densCols(DATA@dr$tsne@cell.embeddings, colramp=colorRampPalette(c("grey80","grey80","navy", "dark green","orange", "red", "firebrick")), nbin = 500),las=1,xlab="tSNE_1",ylab="tSNE_2")
dev.off()
for(k in 1:10){
  png(filename = paste0(opt$output_path,"/tSNE_plots/tSNE_density_",k,".png"),width = 600,height =650,res = 100)
  plot(DATA@dr$tsne@cell.embeddings,pch=16,cex=0.5, col=densCols(DATA@dr$tsne@cell.embeddings, colramp=colorRampPalette(c("grey80","grey80","navy", "dark green","orange", "red", "firebrick")), nbin = 500,bandwidth = k),las=1,xlab="tSNE_1",ylab="tSNE_2")
  dev.off() }
#---------




### Finding clusters using SNN
#---------
cat("\nClustering with SNN ...\n")
if(!dir.exists(paste0(opt$output_path,"/clustering_SNN"))){dir.create(paste0(opt$output_path,"/clustering_SNN"))}
for(k in seq(.05,2,by=.05)){
  if(k!=.05){
    DATA <- FindClusters(object = DATA, reduction.type = "pca", dims.use = 1:5, resolution = k, print.output = F)
  } else { DATA@ident <- factor(NULL) ; DATA <- FindClusters(object = DATA, reduction.type = "pca", dims.use = 1:top_PCs, resolution = k, print.output = F, save.SNN = TRUE,force.recalc = T)}
  png(filename = paste0(opt$output_path,"/clustering_SNN/tSNE_res.",k,".png"),width = 700,height = 600,res = 150)
  TSNEPlot(object = DATA, group.by=paste0("res.",k), pt.size = .5, plot.title= paste0("Clustering (res.",k,")"))
  dev.off()
}
#---------




### Clustering using HDBSCAN
#---------
cat("\nClustering with HDBSCAN ...\n")
if(!dir.exists(paste0(opt$output_path,"/clustering_HDBSCAN"))){dir.create(paste0(opt$output_path,"/clustering_HDBSCAN"))}
for(k in seq(5,100,by=2)){
  clusters <- hdbscan(DATA@dr$tsne@cell.embeddings,minPts = k)
  names(clusters$cluster) <- rownames(DATA@meta.data)
  DATA <- AddMetaData(object = DATA, metadata = clusters$cluster, col.name = paste0("hdbscan.",k))
  png(filename = paste0(opt$output_path,"/clustering_HDBSCAN/tSNE_hdbscan.",k,".png"),width = 700,height = 600,res = 150)
  TSNEPlot(object = DATA, group.by=paste0("hdbscan.",k), pt.size = .5, plot.title= paste0("Clustering (hdbscan.",k,")"))
  dev.off()
}
#---------




### Clustering using FLOWPEAKS
#---------
cat("\nClustering with FLOWPEAKS ...\n")
if(!dir.exists(paste0(opt$output_path,"/clustering_FlowPeaks"))){dir.create(paste0(opt$output_path,"/clustering_FlowPeaks"))}
for(k in seq(.05,4,by=.05)){
  c <- flowPeaks(DATA@dr$tsne@cell.embeddings,h0 = k)
  cluster_ID <- assign.flowPeaks(c,DATA@dr$tsne@cell.embeddings,tol = 0.01,fc = .1)
  names(cluster_ID) <- rownames(DATA@meta.data)
  DATA <- AddMetaData(object = DATA, metadata = cluster_ID, col.name = paste0("flowpeaks.",k))
  png(filename = paste0(opt$output_path,"/clustering_FlowPeaks/tSNE_flowpeaks.",k,".png"),width = 700,height = 600,res = 150)
  TSNEPlot(object = DATA, group.by=paste0("flowpeaks.",k), pt.size = .3, plot.title= paste0("Clustering (flowpeaks.",k,")"))
  dev.off()
}
#---------




### Clustering using DBSCAN
#---------
cat("\nClustering with DBSCAN ...\n")
if(!dir.exists(paste0(opt$output_path,"/clustering_DBSCAN"))){dir.create(paste0(opt$output_path,"/clustering_DBSCAN"))}
e <- c()
for(k in seq(2,30,by=.1)){
  clusters <- dbscan(DATA@dr$tsne@cell.embeddings, eps = k ,minPts = 30)
  e <- c(e,table(clusters$cluster)[1]/length(clusters$cluster))
  names(clusters$cluster) <- rownames(DATA@meta.data)
  DATA <- AddMetaData(object = DATA, metadata = clusters$cluster, col.name = paste0("dbscan.",k))
  png(filename = paste0(opt$output_path,"/clustering_DBSCAN/tSNE_dbscan.",k,".png"),width = 700,height = 600,res = 150)
  TSNEPlot(object = DATA, group.by=paste0("dbscan.",k), pt.size = .3, plot.title= paste0("Clustering (dbscan.",k,")"))
  dev.off()
}
#---------



### Saving the Seurat object
#---------
write.csv(DATA@meta.data,paste0(opt$output_path,"/Metadata_with_clustering.csv"),row.names = T)
saveRDS(DATA, file = paste0(opt$output_path,"/Seurat_object.rds") )
#---------


cat("\n!!! Script executed Sucessfully !!!\n")


### System and session information
#---------
cat("\n\n\n\n... SYSTEM INFORMATION ...\n")
Sys.info()

cat("\n\n\n\n... SESSION INFORMATION ...\n")
sessionInfo()
#---------