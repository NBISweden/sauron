#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
cat("\nRunning DIFFERENTIAL EXPRESSION with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object FILE."),
  make_option(c("-c", "--clustering_use"),        type = "character",   metavar="character",   default='none',  help="The clustering column to be used. Should be chosen from one of the method in script 02."),
  make_option(c("-m", "--metadata_use"),          type = "character",   metavar="character",   default='none',  help="A metadata column name can be supplied to be compared for differential expression for each clusters, as a second step. If nothing or incorrect values are provided, this step is ignored. Multible parameters can be provided and should be comma separated: 'TERM1,TERM2'. In this case, the analysis will be done for each factor individually."),
  make_option(c("-e", "--exclude_cluster"),       type = "character",   metavar="character",   default='none',  help="Clusters to be exluded from the analysis. Usefull for removing outlier cells. If the cluster name is not found, this will be ignored. Multible parameters can be provided and should be comma separated: 'TERM1,TERM2'. In this case, both clusters will be excluded from the DGE analysis."),
  make_option(c("-n", "--DGE_method"),            type = "character",   metavar="character",   default='wilcox',  help='Method to be used for differential gene expression (as implemented in Seurat). Available options are: "wilcox" (default),"bimod","roc","t","MAST","LR","DESeq2","negbinom","poisson".'),
  make_option(c("-p", "--only_positive"),         type = "character",   metavar="character",   default='true',  help='Whether to return only positive fold changes'),
  make_option(c("-l", "--covariates"),            type = "character",   metavar="character",   default='NULL',  help='Which covariates to include in the differential expression (defaul is none). Should be names in the metadata. Mutiple can be used comma separated.'),
  make_option(c("-d", "--max_cells_per_ident"),   type = "character",   metavar="character",   default='100',  help='The maximum number of cells to be included in the comparisson. Cells are sampled randomly.'),
  make_option(c("-f", "--pvalue_threshold"),      type = "character",   metavar="character",   default='0.1',  help='The minimum value for pvalue to return.'),
  make_option(c("-j", "--logfc_threshold"),       type = "character",   metavar="character",   default='0.1',  help='The minimum value for the absolute logFC to return.'),
  make_option(c("-a", "--assay"),                 type = "character",   metavar="character",   default='RNA',  help="Assay to be used in the analysis."),
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
pkgs <- c("rafalib","dplyr","RColorBrewer","scales","igraph","cowplot","ggplot2","Seurat")

library(Seurat)
library(dplyr)
library(scales)
library(RColorBrewer)
library(igraph)
library(rafalib)
library(parallel)
library(cowplot)
library(ggplot2)

#inst_packages(pkgs)
#---------



#############################
### LOAD Seurat.v3 OBJECT ###
#############################
cat("\n### LOADING Seurat.v3 OBJECT ###\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------


###################################################################
### Finding differentially expressed genes (cluster biomarkers) ###
###################################################################
DATA@active.ident <- factor(NULL)
DATA <- SetIdent(DATA,value = factor(as.character(DATA@meta.data[,opt$clustering_use])))

#If the cluster to be excluded is present in the data, it will be removed
if(sum(as.character(unlist(strsplit(opt$exclude_cluster,","))) %in% unique(DATA@meta.data[,opt$clustering_use])) > 0 ){
  DATA <- SubsetData(DATA, cells = colnames(DATA)[! (DATA@meta.data[,opt$clustering_use] %in% as.character(unlist(strsplit(opt$exclude_cluster,","))) )]) #Filter out cell with no assigned clusters
  DATA@meta.data <- DATA@meta.data[colnames(DATA)[! (DATA@meta.data[,opt$clustering_use] %in% as.character(unlist(strsplit(opt$exclude_cluster,","))) )],]
}

if(file.exists(paste0(opt$output_path,"/Cluster_marker_genes.csv"))){
  cat("\nThe following Marker gene list was found here:\n",paste0(opt$output_path,"/Cluster_marker_genes.csv"),"\nDifferential expression among clusters will be skiped\n")
  DATA_markers <- droplevels(read.csv2(paste0(opt$output_path,"/Cluster_marker_genes.csv"),row.names = 1))
  
} else {
  cat("\nThe Marker gene list was not found in the output folder. Computing differential expression among clusters ...")

  cat("\nThe following clusters (",opt$clustering_use,") will be used for analysis ...\n")
  print(as.character(unique(DATA@meta.data[,opt$clustering_use])))
  
  #DATA <- BuildSNN(DATA,reduction.type = "tsne",plot.SNN = F,k.param = 3,prune.SNN = .1)
  DATA_markers <- FindAllMarkers(object = DATA, assay = opt$assay, only.pos = as.logical(opt$only_positive),
                                 min.pct = 0.1, min.diff.pct = 0.05,max.cells.per.ident = as.numeric(opt$max_cells_per_ident),print.bar = T,
                                 do.print = T,return.thresh = as.numeric(opt$pvalue_threshold),test.use = opt$DGE_method,
                                 logfc.threshold = as.numeric(opt$logfc_threshold))

  write.csv2(DATA_markers,file = paste0(opt$output_path,"/Cluster_marker_genes.csv"),row.names = T)
}
#---------



##########################################################
### Plot a heatmap with the top genes for each cluster ###
##########################################################
cat("\nPlotting heatmap of Cluster Marker genes ...\n")
DATA_markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
#pdf(paste0(opt$output_path,"/Cluster_markers_heatmap.pdf"),width = 10,height = 10, useDingbats = F)
png(filename = paste0(opt$output_path,"/Cluster_markers_heatmap.png"),width = 1500,height = 1500,res = 150)
try(DoHeatmap(object = DATA, features = as.character(unique(top10$gene)), assay=opt$assay))
invisible(dev.off())
#---------



#############################################################
### Plot violin plots with the top genes for each cluster ###
#############################################################
png(filename = paste0(opt$output_path,"/violinPlot_genes_per_cluster_scaled.png"),width = 200*30,height = 200*1.5*length(as.character(unique(top10$gene)))/10,res = 150)
VlnPlot(object = DATA, features = as.character(unique(top10$gene)), pt.size = .1, ncol=10,assay = opt$assay)
invisible(dev.off())

png(filename = paste0(opt$output_path,"/UMAP_plot_genes_per_cluster_scaled.png"),width = 350*30,height = 300*3*length(as.character(unique(top10$gene)))/10,res = 150)
FeaturePlot(object = DATA, features = as.character(unique(top10$gene)), cols = c("grey", "blue"), reduction = "umap",pt.size = 1,ncol=10)
invisible(dev.off())
#---------



if(sum(as.character(unlist(strsplit(opt$metadata_use,","))) %in% colnames(DATA@meta.data)) > 0 ){
  for(k in as.character(unlist(strsplit(opt$metadata_use,",")))){
    
######################################################################
### Plot the percentage of metadata in each cluster and vice-versa ###
######################################################################
out <- paste0(opt$output_path,"/Analysis_",k)
if(!dir.exists(out)){dir.create(out,recursive = T)}
cat("\nPlotting cell distribution across\t",k," ...\n")
my_meta_levels <- as.character(unique(DATA@meta.data[,k]))

proportion <- as.data.frame(lapply(levels(DATA@active.ident),function(x){c(unname(table(DATA@meta.data[DATA@active.ident==x,k]))) [1:length(my_meta_levels)]} ))
proportion[is.na(proportion)] <- 0
rownames(proportion) <- my_meta_levels
colnames(proportion) <- levels(DATA@active.ident)

sa <- cbind(stack(as.data.frame(proportion/rowSums(proportion))), rep(rownames(proportion),ncol(proportion)) )
colnames(sa) <- c("prop","clusters",k)

proportion2 <- t(as.matrix(t(proportion)/colSums(proportion)))
cl_order <- order(proportion2[my_meta_levels[1],],proportion2[my_meta_levels[2],],decreasing=T)
sa2 <- cbind(stack(as.data.frame(t(t(proportion[,cl_order])/colSums(proportion[,cl_order])) )), rep(rownames(proportion),ncol(proportion)) )
colnames(sa2) <- c("prop","clusters",k)

assign(k,k)

png(filename = paste0(out,"/barplot_overall.png"),width = 1300,height = 500,res = 150)
print(plot_grid(ggplot(data=sa2,mapping = aes(x=clusters, y=prop, fill=get(k)) ) + geom_bar(stat="identity"),ncol = 2,
          ggplot(data=sa,mapping = aes(x=get(k), y=prop, fill=clusters) ) + geom_bar(stat="identity") ))
invisible(dev.off())
#--------- 



####################################################
### Identifying relevant markers across clusters ###
####################################################
cat("Calculating Differential gene expression among ",k,"   for each cluster ...\n")
marker_list <- list()
cluster_data <- list()
for(i in unique(DATA@active.ident)){
  cat(paste0("... Processing cell cluster #",i," ..."),"\n")
  temp <- SubsetData(DATA, cells = colnames(DATA)[DATA@active.ident == i]) #Select cell from a cluster
  temp@active.ident <- factor(NULL)
  temp <- SetIdent(temp,value = temp@meta.data[,k])
  
  #check if all conditions have cells for differential expression
  if( length(table(temp@active.ident)) == length(my_meta_levels)  ){
    
    #check if the differential expression was previously calculated and saved as a .csv file
    if(file.exists(paste0(out,"/DEGs_in_cluster",i,".csv"))){
      
      cat("The following DEG list was found here:\n",paste0(out,"/DEGs_in_cluster",i,".csv"),"\nDifferential expression among clusters will be skiped\n")
      temp_markers <- data.frame(read.csv2(paste0(out,"/DEGs_in_cluster",i,".csv"),row.names = 1))
      
      } else {
        cat("\nThe following DEG list was not found in the output folder. Computing differential expression ...")
      
        try( temp_markers <- FindAllMarkers(object = temp, only.pos = T, assay = opt$assay) )
        try( temp_markers <- temp_markers[(temp_markers$p_val < 0.01)&(temp_markers$avg_logFC > 0.3),] )
        try( temp_markers <- temp_markers[order(temp_markers$p_val),] )
        try( marker_list[[i]] <- temp_markers )
        try( cluster_data[[i]] <- as.matrix(temp@assays[[opt$assay]]@data)[temp_markers$gene,] )
        try( write.csv2(temp_markers,paste0(out,"/DEGs_in_cluster",i,".csv"),row.names = T) )
        #if(nrow(temp_markers) > 20){top_temp <- temp_markers[1:20,]}
        
        }
      
      try( temp_markers %>% group_by(cluster) %>% top_n(15, avg_logFC) -> top_temp)
      
      #Print the top differentially expressed genes
      png(filename = paste0(out,"/DEGs_in_cluster",i,".png"),width = 200*10*2,height = 200*1.5*length(as.character(unique(top_temp$gene)))/10,res = 150)
      print(VlnPlot(object = temp, features = as.character(unique(top_temp$gene)), pt.size = .1,ncol = 10,assay = opt$assay)) #it does not work if you don't have the print command in front of it!
      dev.off()
      
    } else {
      cat(paste0("... Cell cluster #",i," has no cells in at least one of the conditions ..."),"\n")
    }
}}}
#---------    



#############################
### SYSTEM & SESSION INFO ###
#############################
#---------
print_session_info()
#---------
