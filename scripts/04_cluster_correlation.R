#!/usr/bin/env Rscript

### LOAD OPTPARSE
#---------
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')};
library(optparse)
#---------



### DEFINE PATH TO LOCAL FILES
#---------
cat("\nRunning MARKER IDENTIFICATION and DIFFERENTIAL EXPRESSION with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object FILE."),
  make_option(c("-c", "--clustering_use"),        type = "character",   metavar="character",   default='none',  help="The clustering column to be used. Should be chosen from one of the method in script 02."),
  make_option(c("-e", "--exclude_cluster"),       type = "character",   metavar="character",   default='none',  help="Clusters to be exluded from the analysis. Usefull for removing outlier cells. If the cluster name is not found, this will be ignored. Multible parameters can be provided and should be comma separated: 'TERM1,TERM2'. In this case, both clusters will be excluded from the DGE analysis."),
  make_option(c("-m", "--merge_cluster"),       type = "character",   metavar="character",   default='0.9,0.8,0.7',  help="Correlation threshold between clusters to be used to merge clusters."),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output DIRECTORY.")
) 
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)
#---------



### LOAD LIBRARIES
#---------
cat("\nLoading/installing libraries ...\n")
initial.options <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(sub("--file=","",initial.options[grep("--file=",initial.options)]))
source( paste0(script_path,"/inst_packages.R") )
pkgs <- c("rafalib","dplyr","RColorBrewer","scales","igraph","pheatmap","Seurat","fields","data.table")
inst_packages(pkgs)
#---------




### LOAD Seurat OBJECT 
#---------
cat("\nLoading/ data and metadata ...\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------



### Finding differentially expressed genes (cluster biomarkers)
#---------
DATA@ident <- factor(NULL)
DATA <- SetIdent(DATA,ident.use = DATA@meta.data[,opt$clustering_use])

#If the cluster to be excluded is present in the data, it will be removed
if(sum(as.character(unlist(strsplit(opt$exclude_cluster,","))) %in% unique(DATA@meta.data[,opt$clustering_use])) > 0 ){
  DATA <- SubsetData(DATA, cells.use = DATA@cell.names[! (DATA@meta.data[,opt$clustering_use] %in% as.character(unlist(strsplit(opt$exclude_cluster,","))) )]) #Filter out cell with no assigned clusters
  DATA@meta.data <- DATA@meta.data[DATA@cell.names[! (DATA@meta.data[,opt$clustering_use] %in% as.character(unlist(strsplit(opt$exclude_cluster,","))) )],]
}
#---------




###Correlation between cluster mean expression
#---------
cat("\nCopputing cluster average correlation ...\n")
avg <- t(rowsum(t(as.matrix(DATA@data)) , group = DATA@ident))
avg <- avg[ rowSums(avg  > .5 ) > 1, ]
dim(avg)
avg <- t(t(avg) / as.numeric(table(DATA@ident)) )
hvgs <- DATA@var.genes[ DATA@var.genes %in% rownames(avg)]
cors <- cor(avg[hvgs,],method = "pearson")

cat("\nPlotting ...\n")
png(filename = paste0(opt$output_path,"/Average_Cluster_correlation_heatmap.png"),width = 700,height = 600,res = 150)
pheatmap(cors,cluster_rows = F,cluster_cols = F,colorRampPalette(c("grey95","grey80","firebrick"))(50),scale = "none",border="white",display_numbers = T)
invisible(dev.off())
#---------




###Correlation between cluster averages to every cell
#---------
cat("\nCopputing per-cell to cluster correlation ...\n")
sc_data <- as.matrix(DATA@data)
hvgs <- DATA@var.genes[ DATA@var.genes %in% rownames(avg)]

#Compute correlations
cor_data <- list()
for( i in sort(as.character(unique(DATA@ident))) ){
  temp <- apply(sc_data,2,function(x){ cor(avg[hvgs,i], x[hvgs],method = "pearson",use = "complete.obs") })
  cor_data[[i]] <- temp
}

#Plotting
cat("\nPlotting ...\n")
n <- length(unique(as.character(DATA@ident)))
png(filename = paste0(opt$output_path,"/Single_cell_Cluster_correlation_heatmap.png"),width = 400*4,height = 350*ceiling(n / 4),res = 150)
par(mar=c(1.5,1.5,3,5), mfrow=c(ceiling(n / 4),4))
for( i in sort(unique(as.character(DATA@ident))) ){
  myCorGrad <- colorRampPalette(c("gray85","gray85","gray70","orange3","firebrick","red"))(10)
  lim <- max(cor_data[[i]]^2)
  temp_cor_data1 <- ((cor_data[[i]]^2) - 0) / (lim - 0)
  temp_cor_data1[temp_cor_data1>1] <- 1
  plot(DATA@dr$tsne@cell.embeddings,pch=20,cex=0.8, line=0.5, col=myCorGrad[ round(temp_cor_data1*9)+1], yaxt="n",xaxt="n",xlab="tSNE1",ylab="tSNE2",lwd=0.25, main=paste0("Cor. to Cluster ",i))
  image.plot(1,1,cbind(0,lim),legend.only = T,col = myCorGrad)
}
dev.off()
#---------




###Correlation between cluster averages to every cell
#---------
cat("\nMerging clusters ...\n")
merge_par <- as.numeric(unlist(strsplit(opt$merge,",")))


for(j in merge_par){
  if( j > min(cors) ){
    tcors <- (cors > j)*1
    cell_clust <- DATA@ident
    clust <- rownames(tcors)

    for( i in clust){
      sel <- rownames(tcors)[ tcors[i,] > 0 ]
      cell_clust[cell_clust %in% sel] <- sel[1]
    }
    DATA <- AddMetaData(object = DATA, metadata = cell_clust, col.name = paste0("merged.",j))
    
    png(filename = paste0(opt$output_path,"/tSNE_merged.",j,".png"),width = 700,height = 600,res = 150)
    TSNEPlot(object = DATA, group.by=paste0("merged.",j), pt.size = .5, plot.title= paste0("Clustering (merged.",j,")"))
    dev.off()
  }
}
#---------




### Saving the Seurat object
#---------
saveRDS(DATA, file = opt$Seurat_object_path )
#---------





cat("\n!!! Script executed Sucessfully !!!\n")


### System and session information
#---------
cat("\n\n\n\n... SYSTEM INFORMATION ...\n")
Sys.info()

cat("\n\n\n\n... SESSION INFORMATION ...\n")
sessionInfo()
#---------