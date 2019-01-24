#!/usr/bin/env Rscript

### LOAD LIBRARIES
#---------
library(optparse)
#---------


### DEFINE PATH TO LOCAL FILES
#---------
cat("\nRunning scQC with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object"),
  make_option(c("-c", "--columns_metadata"),      type = "character",   metavar="character",   default='none',  help="Column names in the Metadata matrix (only factors allowed, not continuous variables)"),
  make_option(c("-r", "--regress"),               type = "character",   metavar="character",   default='none',  help="Variables to regress out"),
  make_option(c("-p", "--PCs_use"),               type = "character",   metavar="character",   default='top,5', help="Method and threshold level for selection of significant principal components. The method should be separated from the threshold via a comma. 'top,5' will use the top 5 PCs, which is the default. 'var,1' will use all PCs with variance above 1%."),
  make_option(c("-v", "--var_genes"),             type = "character",   metavar="character",   default='yes',  help="Whether use ('yes') or not ('no') only the variable genes for PCA and tSNE. Defult is 'Yes'."),
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
pkgs <- c("Seurat","rafalib","scran","biomaRt","scater","dplyr","RColorBrewer","dbscan","flowPeaks","scales","igraph")
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




### Filter for genes with lowly-expressed genes
#---------
cat("\nNormalizing and regressing uninteresting factors ...\n")
sel <- rowSums(as.matrix(DATA@raw.data) >= 2) >= 5
DATA <- CreateSeuratObject(as.matrix(DATA@raw.data[sel,cells_use]), meta.data = DATA@meta.data[cells_use,])

cat("\n")
print(DATA)
cat("\n")


DATA <- NormalizeData(DATA)
dup <- sum(duplicated.array(DATA@data, MARGIN = 2))
if( dup > 0 ){
  cat("\n ",dup," duplicated values were found in your data. A very small random floating value (sd=0.0001) will be added to each number to overcome the issue ...\n")
  DATA@data[DATA@data != 0] <- DATA@data[DATA@data != 0] + rnorm(sum(DATA@data != 0),sd = 0.0001)
}


if(opt$var_genes == "no"){
  DATA@var.genes <- rownames(DATA@data)
} else {
  #Defining the variable genes based on the mean gene expression abothe the 5% quantile and the dispersion above 2.
  DATA <- FindVariableGenes(object = DATA, mean.function = ExpMean, dispersion.function = LogVMR, y.cutoff = 2,num.bin = 200)
  m <- max(quantile(DATA@hvg.info$gene.mean,probs = c(.025)) , 0.01)
  DATA <- FindVariableGenes(object = DATA, mean.function = ExpMean, dispersion.function = LogVMR, y.cutoff = 2,num.bin = 200,x.low.cutoff = m)
  write.csv2(DATA@hvg.info, paste0(opt$output_path,"/HVG_info.csv"))
}


#Plotting HVGs
png(filename = paste0(opt$output_path,"/Var_gene_selection.png"),width = 700,height = 750,res = 150)
plot(log2(DATA@hvg.info$gene.mean),DATA@hvg.info$gene.dispersion.scaled,cex=.1,main="HVG selection",
     col=ifelse(rownames(DATA@hvg.info)%in% DATA@var.genes,"red","black" ),ylab="scaled.dispersion",xlab="log2(avg. expression)")
abline(v=log2(m),h=2,lty=2,col="grey20",lwd=1)
dev.off()

cat(as.character(unlist(strsplit(opt$regress,","))),"\n")
DATA <- ScaleData(DATA,vars.to.regress = as.character(unlist(strsplit(opt$regress,","))))
#---------



### Running Dimentionality reduction PCA and tSNE
#---------
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
points( var_expl,bg=ifelse( var_expl > .01, "orange","grey"),pch=21)
invisible(dev.off())


if(file.exists(paste0(opt$output_path,"/tSNE_plots/tSNE_coordinates.csv"))){
  cat("\nPre-computed tSNE found and will be used:\n",paste0(opt$output_path,"/tSNE_plots/tSNE_coordinates.csv"),"\n")
  DATA <- SetDimReduction(DATA, reduction.type = "tsne",  slot = "cell.embeddings", new.data = as.matrix(read.csv2(paste0(opt$output_path,"/tSNE_plots/tSNE_coordinates.csv"),row.names = 1)) )
  DATA <- SetDimReduction(DATA, reduction.type = "tsne",  slot = "key", new.data = "tSNE_")
} else {
  cat("\nPre-computed tSNE NOT found. Computing tSNE ...\n")
  DATA <- RunTSNE(object = DATA, perplexity=30, max_iter=2000,theta=0,eta=2000,exaggeration_factor=12,dims.use = 1:5,verbose = T)
  write.csv2(DATA@dr$tsne@cell.embeddings, paste0(opt$output_path,"/tSNE_plots/tSNE_coordinates.csv"))
}
#---------




#Plotting tSNE plots for metadata
#---------
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
for(k in seq(.05,4,by=.05)){
  if(k!=.05){
    DATA <- FindClusters(object = DATA, reduction.type = "pca", dims.use = 1:5, resolution = k, print.output = F)
  } else { DATA@ident <- factor(NULL) ; DATA <- FindClusters(object = DATA, reduction.type = "pca", dims.use = 1:5, resolution = k, print.output = F, save.SNN = TRUE,force.recalc = T)}
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