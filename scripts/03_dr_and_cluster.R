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
  make_option(c("-p", "--PCs_use"),               type = "character",   metavar="character",   default='top,20', help="Method and threshold level for selection of significant principal components. The method should be separated from the threshold via a comma. 'top,5' will use the top 5 PCs, which is the default. 'var,1' will use all PCs with variance above 1%."),
  make_option(c("-v", "--var_genes"),             type = "character",   metavar="character",   default='Seurat,1.5',  help="Whether use 'Seurat' or the 'Scran' method for variable genes identification. An additional value can be placed after a comma to define the level of dispersion wanted for variable gene selection. 'Seurat,2' will use the threshold 2 for gene dispersions. Defult is 'Seurat,1.5'. For Scran, the user should inpup the level of biological variance 'Scran,0.2'. An additional blocking parameter (a column from the metadata) can ba supplied to 'Scran' method block variation comming from uninteresting factors, which can be parsed as 'Scran,0.2,Batch'."),
  make_option(c("-s", "--cluster_use"),           type = "character",   metavar="character",   default='all',    help="The cluster of cells to be used for analysis. Should be defined as the clustering name followed by the cluster names to be used, comma-separated. E.g.: 'SNN_0.2,1,2,3,5,6'."),
  make_option(c("-m", "--cluster_method"),        type = "character",   metavar="character",   default='snn,hc', help="The clustering method and cluster to select for analysis. Current methods are 'snn','dbscan','hdbscan','flowpeaks'. If no input is suplied, all methods will be run."),
  make_option(c("-d", "--dim_reduct_use"),        type = "character",   metavar="character",   default='umap',  help="Which dimensionality reduction method to be run on top of PCA: UMAP (default) or tSNE. If both, then specify them comma-separated'UMAP,tSNE'."),
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
if(!dir.exists(paste0(opt$output_path,"/variable_genes"))){dir.create(paste0(opt$output_path,"/variable_genes"),recursive = T)}
#---------



##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")
initial.options <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(sub("--file=","",initial.options[grep("--file=",initial.options)]))
source( paste0(script_path,"/inst_packages.R") )
source( paste0(script_path,"/compute_hvgs.R") )
pkgs <- c("Seurat","rafalib","scran","biomaRt","scater","dplyr","RColorBrewer","dbscan","flowPeaks","scales","igraph","sva")
inst_packages(pkgs)
#---------



#############################
### LOAD Seurat.v3 OBJECT ###
#############################
cat("\n### LOADING Seurat.v3 OBJECT ###\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------



###################################
### SELECT CELLS FROM A CLUSTER ###
###################################
cat("\n### SELECTING CELLS FROM A CLUSTER ###")
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
    cat("\nThe name of the cluster or the cluster name were not found in your data.\n All cells will be used ...\n")
    cells_use <- rownames(DATA@meta.data)
  }
sel <- rowSums(as.matrix(DATA@assays$RNA@counts) >= 1) >= 1
DATA <- CreateSeuratObject(as.matrix(DATA@assays$RNA@counts[sel,cells_use]), meta.data = DATA@meta.data[cells_use,])
#---------



###########################
### FIND VARIABLE GENES ###
###########################
output_path <- paste0(opt$output_path,"/variable_genes")
if(!dir.exists(output_path)){dir.create(output_path,recursive = T)}
DATA <- compute_hvgs(DATA,VAR_choice,output_path)
#---------



#############################################
### Scaling data and regressing variables ###
#############################################
cat("\n### Scaling data and regressing uninteresting factors ###\n")
vars <- as.character(unlist(strsplit(opt$regress,",")))
DATA <- ScaleData(DATA,vars.to.regress = vars, assay = DefaultAssay(DATA))
#---------

# a@reductions$umap@cell.embeddings
# a <- CreateDimReducObject(embeddings = a@reductions$umap@cell.embeddings,key = "lll")
# 
# a <- compute_hvgs(DATA,VAR_choice,output_path)
# a <- ScaleData(DATA,features = paste0("dim",1:myinput$d), assay = DefaultAssay(DATA))
# a <- RunPCA(a,features = paste0("dim",1:myinput$d))
# a <- RunUMAP(a,dim=1:20)
# UMAPPlot(a,group.by="orig.ident")

###################
### Running PCA ###
###################
cat("\n### Running PCA ###\n")
if(!dir.exists(paste0(opt$output_path,"/PCA_plots"))){dir.create(paste0(opt$output_path,"/PCA_plots"),recursive = T)}
DATA <- RunPCA(DATA, do.print = F, pcs.compute = 50, assay = DefaultAssay(DATA))

ggsave(PCHeatmap(DATA,ncol=5,dims=1:10),filename = paste0("PCA_heatmap.png"), path = paste0(opt$output_path,"/PCA_plots"), dpi = 300,units = "mm",width = 150*5,height = 150*4.5)

var_expl <- (DATA@reductions$pca@stdev^2)/sum(DATA@reductions$pca@stdev^2)
PC_choice <- as.character(unlist(strsplit(opt$PCs_use,",")))
if(casefold(PC_choice[1]) == "var"){  top_PCs <- sum( var_expl > as.numeric(PC_choice[2])/100 )
} else if(casefold(PC_choice[1]) == "top"){top_PCs <- as.numeric(PC_choice[2])} 

png(filename = paste0(opt$output_path,"/PCA_plots/Varianc_explained_PC.png"),width = 1500,height =1200,res = 200)
plot( var_expl,yaxs="i",bg="grey",pch=21,type="l",ylab="% Variance",xlab="PCs",main="all cells",las=1)
points( var_expl,bg=c(rep("orange",top_PCs),rep("grey",50-top_PCs)),pch=21)
invisible(dev.off())
#---------



####################
### Running tSNE ###
####################
if( "tsne" %in% casefold(unlist(strsplit(opt$dim_reduct_use,",")))){
  cat("\n### Running BH-tSNE ###\n")
  if(!dir.exists(paste0(opt$output_path,"/tSNE_plots"))){dir.create(paste0(opt$output_path,"/tSNE_plots"),recursive = T)}
  
  if(file.exists(paste0(opt$output_path,"/tSNE_plots/tSNE_coordinates.csv"))){
    cat("\nPre-computed tSNE found and will be used:\n",paste0(opt$output_path,"/tSNE_plots/tSNE_coordinates.csv"),"\n")
    DATA@reductions[["tsne"]] <- CreateDimReducObject(embeddings = as.matrix(read.csv2(paste0(opt$output_path,"/tSNE_plots/tSNE_coordinates.csv"),row.names = 1)),key = "tSNE_",assay = DefaultAssay(DATA))
    } else { cat("\nPre-computed tSNE NOT found. Computing tSNE ...\n")
    DATA <- RunTSNE(object = DATA, perplexity=30, max_iter=2000,theta=0.1,eta=2000,exaggeration_factor=12,dims.use = 1:top_PCs,verbose = T,num_threads=0)
    write.csv2(DATA@reductions$tsne@cell.embeddings, paste0(opt$output_path,"/tSNE_plots/tSNE_coordinates.csv"))}
}
#---------



####################
### Running UMAP ###
####################
if( "umap" %in% casefold(unlist(strsplit(opt$dim_reduct_use,",")))){
  cat("\n### Running UMAP ###\n")
  if(!dir.exists(paste0(opt$output_path,"/UMAP_plots"))){dir.create(paste0(opt$output_path,"/UMAP_plots"),recursive = T)}
  
  if(file.exists(paste0(opt$output_path,"/UMAP_plots/UMAP_coordinates.csv"))){
    cat("\nPre-computed UMAP found and will be used:\n",paste0(opt$output_path,"/UMAP_plots/UMAP_coordinates.csv"),"\n")
    DATA@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(read.csv2(paste0(opt$output_path,"/UMAP_plots/UMAP_coordinates.csv"),row.names = 1)),key = "UMAP_",assay = DefaultAssay(DATA))
  } else { cat("\nPre-computed UMAO NOT found. Computing UMAP ...\n")
    DATA <- RunUMAP(object = DATA, n.neighbors = 50, dims = 1:top_PCs,min.dist = 0.001, n.components = 2, verbose = T,num_threads=0)
    DATA <- RunUMAP(object = DATA, n.neighbors = 50, dims = 1:top_PCs,min.dist = 0.001, n.components = 10, verbose = T,num_threads=0,reduction.name = "UMAP10",reduction.key = "umap10_")
    write.csv2(DATA@reductions$umap@cell.embeddings, paste0(opt$output_path,"/UMAP_plots/UMAP_coordinates.csv"))}
}
#---------



###################
### Running SNN ###
###################
cat("\n### Running SNN ###\n")
DATA <- FindNeighbors(DATA,assay = DefaultAssay(DATA),graph.name="SNN")
g <- graph_from_adjacency_matrix(as.matrix(DATA@graphs$SNN),weighted = T,diag = F,mode = "undirected")

for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
  png(filename = paste0(opt$output_path,"/",i,"_plots","/",i,"_plot_with_SNN_overlay.png"),width = 200,height =205,res = 600,units = "mm")
  plot.igraph(x = g, layout = DATA@reductions[[i]]@cell.embeddings[,1:2], edge.width = E(graph = g)$weight/4, vertex.label = NA,
              edge.color = colorRampPalette(c("grey90","black"))(50)[round(E(graph = g)$weight*49+1)],
              vertex.size = 1,vertex.frame.color=hue_pal(l=50, c=80)(length(levels(factor(DATA$orig.ident))))[factor(DATA$orig.ident)],
              vertex.color = hue_pal()(length(unique(DATA$orig.ident)))[factor(DATA$orig.ident)])
  invisible(dev.off())
}
#---------



#########################################
### Plotting Dimensionality Reduction ###
#########################################
mtdt <- c("nCount_RNA","nFeature_RNA","S.Score","G2M.Score","percent_rpl","percent_rps","percent_mt-")
mtdt <- mtdt[mtdt %in% colnames(DATA@meta.data)]
j <- as.character(unlist(strsplit(opt$columns_metadata,",")))

for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
  temp <- FeaturePlot(object = DATA, features = mtdt, cols = col_scale,pt.size = .5,reduction = i,ncol = 5,dims = 1:2)
  ggsave(temp,filename = paste0(i,"_metadata_dim1_dim2.png"), path = paste0(opt$output_path,"/",i,"_plots"), dpi = 300,units = "mm",width = 170*5,height = 150*ceiling(length(mtdt)/5) )
  
  temp2 <- DimPlot(DATA,dims = 1:2,reduction = i,group.by = j,pt.size = .3,ncol = 5)
  ggsave(temp2,filename = paste0(i,"_metadata_factors_dim1_dim2.png"), path = paste0(opt$output_path,"/",i,"_plots"), dpi = 300,units = "mm",width = 170*5,height = 150*ceiling(length(j)/5) )
  
  if(i == "pca"){
    temp <- FeaturePlot(object = DATA, features = mtdt, cols = col_scale,pt.size = .5,reduction = i,ncol = 5,dims = 3:4)
    ggsave(temp,filename = paste0(i,"_metadata_dim3_dim4.png"), path = paste0(opt$output_path,"/",i,"_plots"), dpi = 300,units = "mm",width = 170*5,height = 150*ceiling(length(mtdt)/5) )
    
    temp2 <- DimPlot(DATA,dims = 3:4,reduction = i,group.by = j,pt.size = .3,ncol = 5)
    ggsave(temp2,filename = paste0(i,"_metadata_factors_dim3_dim4.png"), path = paste0(opt$output_path,"/",i,"_plots"), dpi = 300,units = "mm",width = 170*5,height = 150*ceiling(length(j)/5) )
} }
#---------



################################
### Clustering using Louvain ###
################################
if( 'snn' %in% casefold(unlist(strsplit(opt$cluster_method,split = ","))) ){
  cat("\n### Clustering with SNN ###\n")
  if(!dir.exists(paste0(opt$output_path,"/clustering"))){dir.create(paste0(opt$output_path,"/clustering"))}
  for(k in seq(.05,2,by=.05)){
    DATA <- FindClusters(object = DATA, reduction.type = "pca", dims.use = 1:top_PCs, resolution = k, verbose = F,graph.name = "SNN")
  }
  for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
  temp2 <- DimPlot(DATA,dims = 1:2,reduction = i,group.by = sort(colnames(DATA@meta.data)[grep("SNN",colnames(DATA@meta.data))]), pt.size = .3,ncol = 8)
  ggplot2::ggsave(temp2,filename = paste0("clustering_SNN_",i,".png"), path = paste0(opt$output_path,"/clustering"), dpi = 300,units = "mm",width = 140*8,height = 100*ceiling(length(grep("SNN",colnames(DATA@meta.data)))/8),limitsize = FALSE )
  rm(temp2); gc()}}
#---------



######################################
### Hierachical clustering on UMAP ###
######################################
if( 'hc' %in% casefold(unlist(strsplit(opt$cluster_method,split = ","))) ){
  cat("\n### Clustering with HC (Hierachical CLustering on UMAP-10dims) ###\n")
  if(!dir.exists(paste0(opt$output_path,"/clustering"))){dir.create(paste0(opt$output_path,"/clustering"))}
  for(k in 2:50){
    cl <- cutree(hclust(dist(DATA@reductions$UMAP10@cell.embeddings,method = "euclidean"),method = "single"),k = k)
    DATA <- AddMetaData(DATA,metadata = setNames(cl,colnames(DATA)), paste0("HC_",k))
  }
  for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
    temp2 <- DimPlot(DATA,dims = 1:2,reduction = i,group.by = colnames(DATA@meta.data)[grep("HC_",colnames(DATA@meta.data))], pt.size = .3,ncol = 8)
    ggplot2::ggsave(temp2,filename = paste0("clustering_HC_",i,".png"), path = paste0(opt$output_path,"/clustering"), dpi = 300,units = "mm",width = 140*8,height = 100*ceiling(length(grep("HC_",colnames(DATA@meta.data)))/8),limitsize = FALSE )
    rm(temp2); gc()}
}
#---------



###################################
### SAVING RAW Seurat.v3 OBJECT ###
###################################
cat("\n### Saving the RAW Seurat object ###\n")
write.csv2(DATA@meta.data,paste0(opt$output_path,"/Metadata_with_clustering.csv"),row.names = T)
saveRDS(DATA, file = paste0(opt$output_path,"/Seurat_object.rds") )
#---------



#############################
### SYSTEM & SESSION INFO ###
#############################
#---------
print_session_info()
#---------


