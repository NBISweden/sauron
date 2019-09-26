#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
cat("\nRunning CLUSTER CORRELATION with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object FILE."),
  make_option(c("-c", "--clustering_use"),        type = "character",   metavar="character",   default='none',  help="The clustering column to be used. Should be chosen from one of the method in script 02."),
  #make_option(c("-a", "--avg_expression_table"),  type = "character",   metavar="character",   default='none',  help="OPTIONAL. Path to a .csv expression matrix contating average expression of cell populations of interest. Each colum can be a sample."),
  make_option(c("-e", "--exclude_cluster"),       type = "character",   metavar="character",   default='none',  help="Clusters to be exluded from the analysis. Usefull for removing outlier cells. If the cluster name is not found, this will be ignored. Multible parameters can be provided and should be comma separated: 'TERM1,TERM2'. In this case, both clusters will be excluded from the DGE analysis."),
  make_option(c("-m", "--merge_cluster"),         type = "character",   metavar="character",   default='0.9,0.8,0.7',  help="Correlation threshold between clusters to be used to merge clusters."),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output DIRECTORY.")
) 
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)
#---------



##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")
initial.options <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(sub("--file=","",initial.options[grep("--file=",initial.options)]))
source( paste0(script_path,"/inst_packages.R") )
pkgs <- c("Seurat","rafalib","biomaRt","scales","fields","pheatmap","ggplot2")

library(Seurat)
library(dplyr)
library(scales)
library(pheatmap)
library(biomaRt)
library(igraph)
library(fields)
library(rafalib)
library(ggplot2)
#inst_packages(pkgs)
#---------



##########################
### LOAD Seurat OBJECT ###
##########################
cat("\nLoading/ data and metadata ...\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------


###################################################################
### Finding differentially expressed genes (cluster biomarkers) ###
###################################################################
DATA@active.ident <- factor(NULL)
DATA <- SetIdent(DATA, value = DATA@meta.data[,opt$clustering_use])

#If the cluster to be excluded is present in the data, it will be removed
if(sum(as.character(unlist(strsplit(opt$exclude_cluster,","))) %in% unique(DATA@meta.data[,opt$clustering_use])) > 0 ){
  DATA <- SubsetData(DATA, cells = colnames(DATA)[! (DATA@meta.data[,opt$clustering_use] %in% as.character(unlist(strsplit(opt$exclude_cluster,","))) )]) #Filter out cell with no assigned clusters
  DATA@meta.data <- DATA@meta.data[colnames(DATA)[! (DATA@meta.data[,opt$clustering_use] %in% as.character(unlist(strsplit(opt$exclude_cluster,","))) )],]
}
#---------



###Correlation between cluster mean expression
#---------
cat("\nCopputing cluster average correlation ...\n")
avg <- t(rowsum(t(as.matrix(DATA@assays[[DATA@active.assay]]@data)) , group = DATA@active.ident))
avg <- avg[ rowSums(avg  > .5 ) > 1, ]
dim(avg)
avg <- t(t(avg) / as.numeric(table(DATA@active.ident)) )
hvgs <- DATA@assays[[DATA@active.assay]]@var.features[ DATA@assays[[DATA@active.assay]]@var.features %in% rownames(avg)]
cors <- cor(avg[hvgs,],method = "pearson")

cat("\nPlotting ...\n")
png(filename = paste0(opt$output_path,"/Average_Cluster_correlation_heatmap.png"),width = 700,height = 600,res = 150)
pheatmap(cors,cluster_rows = F,cluster_cols = F,colorRampPalette(c("grey95","grey80","firebrick"))(50),scale = "none",border="white",display_numbers = T)
invisible(dev.off())
#---------



##########################################################
### Correlation between cluster averages to every cell ###
##########################################################
cat("\nCopputing per-cell to cluster correlation ...\n")
sc_data <- as.matrix(DATA@assays[[DATA@active.assay]]@data)
hvgs <- DATA@assays[[DATA@active.assay]]@var.features[ DATA@assays[[DATA@active.assay]]@var.features %in% rownames(avg)]

#Compute correlations
cor_data <- list()
for( i in sort(as.character(unique(DATA@active.ident))) ){
  temp <- apply(sc_data,2,function(x){ cor(avg[hvgs,i], x[hvgs],method = "pearson",use = "complete.obs") })
  cor_data[[i]] <- temp
}

#Plotting
cat("\nPlotting ...\n")
n <- length(unique(as.character(DATA@active.ident)))
png(filename = paste0(opt$output_path,"/cell_to_cluster_correlation.png"),width = 400*8,height = 350*ceiling(n / 8),res = 150)
par(mar=c(1.5,1.5,3,5), mfrow=c(ceiling(n / 8),8))
for( i in sort(unique(as.character(DATA@active.ident))) ){
  myCorGrad <- colorRampPalette(c("gray85","gray85","gray70","orange3","firebrick","red"))(10)
  lim <- max(cor_data[[i]]^2)
  temp_cor_data1 <- ((cor_data[[i]]^2) - 0) / (lim - 0)
  temp_cor_data1[temp_cor_data1>1] <- 1
  o <- order(temp_cor_data1)
  plot(DATA@reductions$umap@cell.embeddings[o,],pch=20,cex=0.5, line=0.5, col=myCorGrad[ round(temp_cor_data1[o]*9)+1], yaxt="n",xaxt="n",xlab="tSNE1",ylab="tSNE2",lwd=0.25, main=paste0("Cor. to Cluster ",i))
  image.plot(1,1,cbind(0,lim),legend.only = T,col = myCorGrad)
}
invisible(dev.off())
#---------



##########################################
### Merging highly correlated clusters ###
##########################################
cat("\nMerging clusters ...\n")
merge_par <- as.numeric(unlist(strsplit(opt$merge,",")))

for(j in merge_par){
  if( j > min(cors) ){
    tcors <- (cors > j)*1
    cell_clust <- DATA@active.ident
    clust <- rownames(tcors)
    
    for( i in clust){
      sel <- rownames(tcors)[ tcors[i,] > 0 ]
      cell_clust[cell_clust %in% sel] <- sel[1]
    }
    DATA <- AddMetaData(object = DATA, metadata = factor(cell_clust), col.name = paste0("merged_",j))
    
    temp <- UMAPPlot(object = DATA, group.by=paste0("merged_",j), pt.size = .5, plot.title= paste0("Clustering (merged_",j,")"))
    ggsave(temp,filename = paste0("UMAP_merged_",j,".png"), path = opt$output_path, dpi = 300,units = "mm",width = 170,height = 150 )
  }
}
#---------



###############################
### SAVING Seurat.v3 OBJECT ###
###############################
cat("\n### Saving Seurat object ###\n")
saveRDS(DATA, file = opt$Seurat_object_path )
#---------



#############################
### SYSTEM & SESSION INFO ###
#############################
#---------
print_session_info()
#---------
