#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
cat("\nRunnign DIMENSIONLITY REDUCTION AND CLUSTERING with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object"),
  make_option(c("-c", "--columns_metadata"),      type = "character",   metavar="character",   default='none',  help="Column names in the Metadata matrix (only factors allowed, not continuous variables)"),
  make_option(c("-r", "--regress"),               type = "character",   metavar="character",   default='none',  help="Variables to be regressed out using linear modeling."),
  make_option(c("-p", "--PCs_use"),               type = "character",   metavar="character",   default='top,30', help="Method and threshold level for selection of significant principal components. The method should be separated from the threshold via a comma. 'top,5' will use the top 5 PCs, which is the default. 'var,1' will use all PCs with variance above 1%."),
  make_option(c("-v", "--var_genes"),             type = "character",   metavar="character",   default='scran',  help="Whether use 'Seurat' or the 'Scran' method for variable genes identification. An additional value can be placed after a comma to define the level of dispersion wanted for variable gene selection. 'Seurat,2' will use the threshold 2 for gene dispersions. Defult is 'scran'. For Scran, the user should inpup the level of biological variance 'Scran,0.2'. An additional blocking parameter (a column from the metadata) can ba supplied to 'Scran' method block variation comming from uninteresting factors, which can be parsed as 'Scran,0.2,Batch'."),
  make_option(c("-s", "--cluster_use"),           type = "character",   metavar="character",   default='all',    help="The cluster of cells to be used for analysis. Should be defined as the clustering name followed by the cluster names to be used, comma-separated. E.g.: 'louvain_0.2,1,2,3,5,6'."),
  make_option(c("-m", "--cluster_method"),        type = "character",   metavar="character",   default='leiden', help="The clustering method and cluster to select for analysis. Current methods are 'hc','louvain','dbscan','hdbscan','flowpeaks','kmeans','leiden'. If no input is suplied, all methods will be run."),
  make_option(c("-n", "--k_nearest_neighbour"),   type = "character",   metavar="character",   default='30', help="The number of nearest neighbours for graph construction."),
  make_option(c("-d", "--dim_reduct_use"),        type = "character",   metavar="character",   default='umap',  help="Which dimensionality reduction method to be run on top of `pre_dim_reduct`: UMAP (default) or tSNE. If both, then specify them comma-separated'UMAP,tSNE'."),
  make_option(c("-k", "--pre_dim_reduct"),        type = "character",   metavar="character",   default='pca',  help="Which dimensionality reduction method to be run, PCA is run by default. If you run 'MNN', you can use it here."),
  make_option(c("-a", "--assay"),                 type = "character",   metavar="character",   default='RNA',  help="Assay to be used in the analysis."),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output directory")
)
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)

VAR_choice <- as.character(unlist(strsplit(opt$var_genes,",")))
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
source( paste0(script_path,"/fast_ScaleData.R") )
pkgs <- c("Seurat","rafalib","scran","biomaRt","scater","dplyr","RColorBrewer","dbscan","scales","igraph","sva")

suppressMessages(suppressWarnings({
  library(Seurat) 
  library(dplyr)
  library(scales)
  library(RColorBrewer)
  library(biomaRt)
  library(igraph)
  library(sva)
  library(rafalib)
  library(parallel)
  library(scran)
  library(scater)
  library(dbscan)
}))
#inst_packages(pkgs)
#---------



#############################
### LOAD Seurat.v3 OBJECT ###
#############################
cat("\n### LOADING Seurat.v3 OBJECT ###\n")
DATA <- readRDS(opt$Seurat_object_path)
DATA@active.assay <- opt$assay
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
    cells_use <- rownames(DATA@meta.data)[factor(DATA@meta.data[,clustering_use]) %in% clusters_to_select]}   #Filter out cells with no assigned clusters
  DATA <- SubsetData(DATA,assay = opt$assay,cells = cells_use)
} else {
  cat("\nThe name of the cluster or the cluster name were not found in your data.\n All cells will be used ...\n")
  cells_use <- rownames(DATA@meta.data) }
# sel <- rowSums(as.matrix(DATA@assays$RNA@counts) >= 1) >= 1
# DATA <- CreateSeuratObject(as.matrix(DATA@assays$RNA@counts[sel,cells_use]), meta.data = DATA@meta.data[cells_use,])
#---------



###########################
### FIND VARIABLE GENES ###
###########################
if(! (casefold( opt$assay ) %in% c("mnn")) ) {
  cat("\n### Computing variable genes\n")
  output_path <- paste0(opt$output_path,"/variable_genes")
  DATA <- compute_hvgs(DATA,VAR_choice,output_path,assay = opt$assay)
} else { 
  DATA@assays[[opt$assay]]@var.features <- rownames(DATA@assays[[casefold( opt$assay )]]@data)
}
#---------



#############################################
### Scaling data and regressing variables ###
#############################################
cat("\n### Scaling data and regressing uninteresting factors in all assays ###\n")
vars_regress <- unlist(strsplit(opt$regress,","))
for(i in names(DATA@assays)){
  cat("\n### Processing assay: ",i,"###\n")
  if( length(vars_regress[vars_regress%in%colnames(DATA@meta.data)]) == 0 ){ vars_regress <- NULL}
  DATA <- fast_ScaleData(DATA, vars.to.regress = vars_regress, assay=i)
  # DATA <- ScaleData(DATA, vars.to.regress = vars_regress, assay=i)
}
#---------



###################
### Running PCA ###
###################
cat("\n### Running PCA ###\n")
if(!dir.exists(paste0(opt$output_path,"/pca_plots"))){dir.create(paste0(opt$output_path,"/pca_plots"),recursive = T)}
DATA <- RunPCA(DATA, do.print = F, assay = opt$assay,npcs = 100)
write.csv2(DATA@reductions$pca@cell.embeddings, paste0(opt$output_path,"/pca_plots/PCA_coordinates.csv"))
write.csv2(DATA@reductions$pca@feature.loadings, paste0(opt$output_path,"/pca_plots/PCA_feature_loadings.csv"))

ggsave(PCHeatmap(DATA,ncol=5,dims=1:10),filename = paste0("PCA_heatmap.png"), path = paste0(opt$output_path,"/pca_plots"), dpi = 300,units = "mm",width = 150*5,height = 150*4.5,limitsize = FALSE)

var_expl <- (DATA@reductions$pca@stdev^2)/sum(DATA@reductions$pca@stdev^2)
PC_choice <- as.character(unlist(strsplit(opt$PCs_use,",")))
if(casefold(PC_choice[1]) == "var"){  top_PCs <- sum( var_expl > as.numeric(PC_choice[2])/100 )
} else if(casefold(PC_choice[1]) == "top"){top_PCs <- as.numeric(PC_choice[2])}

png(filename = paste0(opt$output_path,"/pca_plots/Varianc_explained_PC.png"),width = 1500,height =1200,res = 200)
plot( var_expl*100,yaxs="i",bg="grey",pch=21,type="l",ylab="% Variance",xlab="PCs",main="all cells",las=1,ylim=c(0,1.2*max(var_expl*100)))
points( var_expl*100,bg=c(rep("orange",top_PCs),rep("grey",100 - top_PCs)),pch=21)
invisible(dev.off())
#---------



####################
### Running tSNE ###
####################
if( "tsne" %in% casefold(unlist(strsplit(opt$dim_reduct_use,",")))){
  cat("\n### Running BH-tSNE ###\n")
  if(!dir.exists(paste0(opt$output_path,"/tsne_plots"))){dir.create(paste0(opt$output_path,"/tsne_plots"),recursive = T)}

  if(file.exists(paste0(opt$output_path,"/tsne_plots/tSNE_coordinates.csv"))){
    cat("\nPre-computed tSNE found and will be used:\n",paste0(opt$output_path,"/tsne_plots/tSNE_coordinates.csv"),"\n")
    DATA@reductions[["tsne"]] <- CreateDimReducObject(embeddings = as.matrix(read.csv2(paste0(opt$output_path,"/tsne_plots/tSNE_coordinates.csv"),row.names = 1)),key = "tSNE_",assay = opt$assay)
  } else { cat("\nPre-computed tSNE NOT found. Computing tSNE ...\n")
    
    ttt <- Sys.time()
    DATA <- RunTSNE(object = DATA, reduction = casefold(opt$pre_dim_reduct), perplexity=30, max_iter=2000,theta=0.1,eta=2000,exaggeration_factor=12,dims.use = 1:top_PCs,verbose = T,num_threads=0)
    cat("multicore tSNE ran in ",difftime(Sys.time(), ttt, units='mins'))
    write.csv2(DATA@reductions$tsne@cell.embeddings, paste0(opt$output_path,"/tsne_plots/tSNE_coordinates.csv"))}
}
#---------



####################
### Running UMAP ###
####################
if( "umap" %in% casefold(unlist(strsplit(opt$dim_reduct_use,",")))){
  cat("\n### Running UMAP on ",casefold(opt$pre_dim_reduct)," ###\n")
  if(!dir.exists(paste0(opt$output_path,"/_plots"))){dir.create(paste0(opt$output_path,"/umap_plots"),recursive = T)}

  if(file.exists(paste0(opt$output_path,"/umap_plots/UMAP_coordinates.csv"))&file.exists(paste0(opt$output_path,"/umap_plots/UMAP10_coordinates.csv"))){
    cat("\nPre-computed UMAP found and will be used:\n",paste0(opt$output_path,"/umap_plots/UMAP_coordinates.csv"),"\n")
    DATA@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(read.csv2(paste0(opt$output_path,"/umap_plots/UMAP_coordinates.csv"),row.names = 1)),key = "UMAP_",assay = opt$assay)
    DATA@reductions[["umap10"]] <- CreateDimReducObject(embeddings = as.matrix(read.csv2(paste0(opt$output_path,"/umap_plots/UMAP10_coordinates.csv"),row.names = 1)),key = "UMAP_",assay = opt$assay)

  } else { cat("\nPre-computed UMAP NOT found. Computing UMAP ...\n")
    
    ttt <- Sys.time()
    DATA <- RunUMAP(object = DATA, reduction = casefold(opt$pre_dim_reduct), dims = 1:top_PCs, n.components = 2, n.neighbors = 10, spread = .3,
                    repulsion.strength = 1, min.dist= .001, verbose = T,num_threads=0,n.epochs = 200,metric = "euclidean",seed.use = 42)
    cat("UMAP_2dimensions ran in ",difftime(Sys.time(), ttt, units='mins'),"\n")
    invisible(gc())
    ttt <- Sys.time()
    
    DATA <- RunUMAP(object = DATA, dims = 1:top_PCs, n.components = 10, n.neighbors = 15, spread = 1, min.dist= .1,metric = "euclidean",verbose = T,num_threads=0,learning.rate = .2,n.epochs = 200,reduction.name = "umap10",reduction.key = "umap10_")
    cat("UMAP_10dimensions ran in ",difftime(Sys.time(), ttt, units='mins'))
    invisible(gc())
    write.csv2(DATA@reductions$umap@cell.embeddings, paste0(opt$output_path,"/umap_plots/UMAP_coordinates.csv"))
    write.csv2(DATA@reductions$umap10@cell.embeddings, paste0(opt$output_path,"/umap_plots/UMAP10_coordinates.csv"))}
  invisible(gc())
}
#---------



#############################
### Running Diffusion Map ###
#############################
if( "dm" %in% casefold(unlist(strsplit(opt$dim_reduct_use,",")))){
  cat("\n### Running Diffusion Map on ",casefold(opt$pre_dim_reduct)," ###\n")
  if(!dir.exists(paste0(opt$output_path,"/dm_plots"))){dir.create(paste0(opt$output_path,"/dm_plots"),recursive = T)}
  
  if(file.exists(paste0(opt$output_path,"/dm_plots/dm_coordinates.csv"))){
    cat("\nPre-computed Diffusion Map found and will be used:\n",paste0(opt$output_path,"/dm_plots/dm_coordinates.csv"),"\n")
    DATA@reductions[["dm"]] <- CreateDimReducObject(embeddings = as.matrix(read.csv2(paste0(opt$output_path,"/dm_plots/dm_coordinates.csv"),row.names = 1)),key = "DC_",assay = opt$assay)
  } else { cat("\nPre-computed Diffusion Map NOT found. Computing Diffusion Map ...\n")
    
    ttt <- Sys.time()
    dm <- destiny::DiffusionMap( DATA@reductions[[casefold(opt$pre_dim_reduct)]]@cell.embeddings[ , 1:top_PCs], k = 20)
    rownames(dm@eigenvectors) <- colnames(DATA)
    DATA@reductions[["dm"]] <- CreateDimReducObject(embeddings = dm@eigenvectors,key = "DC_",assay = opt$assay)
    cat("multicore Diffusion Map ran in ",difftime(Sys.time(), ttt, units='mins'))
    write.csv2(DATA@reductions$dm@cell.embeddings, paste0(opt$output_path,"/dm_plots/dm_coordinates.csv"))}
}
#---------



###################
### Running ICA ###
###################
if( "ica" %in% casefold(unlist(strsplit(opt$dim_reduct_use,",")))){
  cat("\n### Running ICA ###\n")
  if(!dir.exists(paste0(opt$output_path,"/ICA_plots"))){dir.create(paste0(opt$output_path,"/ICA_plots"),recursive = T)}
  
  if(file.exists(paste0(opt$output_path,"/ICA_plots/ICA_coordinates.csv"))){
    cat("\nPre-computed ICA found and will be used:\n",paste0(opt$output_path,"/ICA_plots/ICA_coordinates.csv"),"\n")
    
    DATA@reductions[["ica"]] <- CreateDimReducObject(embeddings = as.matrix(read.csv2(paste0(opt$output_path,"/ICA_plots/ICA_coordinates.csv"),row.names = 1)),key = "ICA_",assay = opt$assay)
  } else { cat("\nPre-computed ICA NOT found. Computing ICA ...\n")
    
    DATA <- RunICA(DATA,assay = opt$assay, nics = 20,reduction.name = "ica")
    write.csv2(DATA@reductions$dm@cell.embeddings, paste0(opt$output_path,"/ICA_plots/ICA_coordinates.csv"))}
}
#---------



###################
### Running SNN ###
###################
cat("\n### Running SNN on ",casefold(opt$pre_dim_reduct)," ###\n")
DATA <- FindNeighbors(DATA, assay = opt$assay, graph.name="SNN", prune.SNN = .2,k.param = as.numeric(opt$k_nearest_neighbour),force.recalc = T,reduction = casefold(opt$pre_dim_reduct), dims = 1:top_PCs )
g <- graph_from_adjacency_matrix(as.matrix(DATA@graphs$SNN),weighted = T,diag=F)
g <- simplify(g)
saveRDS(DATA@graphs$SNN, file = paste0(opt$output_path,"/SNN_Graph.rds") )

for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
  png(filename = paste0(opt$output_path,"/",i,"_plots","/",i,"_plot_with_SNN_overlay.png"),width = 200,height =205,res = 600,units = "mm")
  plot.igraph(x = g, layout = DATA@reductions[[i]]@cell.embeddings[,1:2], edge.width = E(graph = g)$weight/4, vertex.label = NA,
              edge.color = colorRampPalette(c("grey90","black"))(50)[round(E(graph = g)$weight/2*49+1)],edge.arrow.size=E(graph = g)$weight*0,
              vertex.size = 1,vertex.frame.color=hue_pal(l=50, c=80)(length(levels(factor(DATA$orig.ident))))[factor(DATA$orig.ident)],
              vertex.color = hue_pal()(length(unique(DATA$orig.ident)))[factor(DATA$orig.ident)])
  invisible(dev.off())
}
rm(g); invisible(gc())
#---------



#########################################
### Plotting Dimensionality Reduction ###
#########################################
mtdt <- colnames(DATA@meta.data) [ grepl("nFeature|nCount|_index|[.]Score",colnames(DATA@meta.data) ) ]
mtdt <- c(mtdt, "perc_mito" ,"perc_rps","perc_rpl","perc_hb", "perc_protein_coding" ,"perc_lincRNA","perc_snRNA","perc_miRNA","perc_processed_pseudogene",
"perc_unknown","perc_Chr_1","perc_Chr_X","perc_Chr_Y","perc_Chr_MT")
mtdt <- mtdt[mtdt %in% colnames(DATA@meta.data)]
j <- unlist(strsplit(opt$columns_metadata,","))

for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
  temp <- FeaturePlot(object = DATA, features = mtdt, cols = col_scale,pt.size = .5,reduction = i,ncol = 5,dims = 1:2)
  ggsave(temp,filename = paste0(i,"_metadata_dim1_dim2.png"), path = paste0(opt$output_path,"/",i,"_plots"), dpi = 300,units = "mm",width = 170*5,height = 150*ceiling(length(mtdt)/5),limitsize = FALSE )

  temp2 <- DimPlot(DATA,dims = 1:2,reduction = i,group.by = j,pt.size = .3,ncol = 5)
  ggsave(temp2,filename = paste0(i,"_metadata_factors_dim1_dim2.png"), path = paste0(opt$output_path,"/",i,"_plots"), dpi = 300,units = "mm",width = 170*5,height = 150*ceiling(length(j)/5),limitsize = FALSE )

  if(i == "pca"){
    temp <- FeaturePlot(object = DATA, features = mtdt, cols = col_scale,pt.size = .5,reduction = i,ncol = 5,dims = 3:4)
    ggsave(temp,filename = paste0(i,"_metadata_dim3_dim4.png"), path = paste0(opt$output_path,"/",i,"_plots"), dpi = 300,units = "mm",width = 170*5,height = 150*ceiling(length(mtdt)/5),limitsize = FALSE )

    temp2 <- DimPlot(DATA,dims = 3:4,reduction = i,group.by = j,pt.size = .3,ncol = 5)
    ggsave(temp2,filename = paste0(i,"_metadata_factors_dim3_dim4.png"), path = paste0(opt$output_path,"/",i,"_plots"), dpi = 300,units = "mm",width = 170*5,height = 150*ceiling(length(j)/5),limitsize = FALSE )
  } }
rm(temp,temp2); invisible(gc())
#---------



################################
### Clustering using Louvain ###
################################
if(  'louvain' %in% casefold(unlist(strsplit(opt$cluster_method,split = ",")))  ){
  cat("\n### Clustering with louvain ###\n")
  if(!dir.exists(paste0(opt$output_path,"/clustering"))){dir.create(paste0(opt$output_path,"/clustering"))}
  #remove old clustering with this name
  DATA@meta.data <- DATA@meta.data[ , ! grepl( "louvain_", colnames(DATA@meta.data) ) ]
  
  DATA <- FindClusters(object = DATA, reduction.type = "pca", dims.use = 1:top_PCs, resolution = seq(.05,2,by=.05), verbose = T,graph.name = "SNN",algorithm = 1)
  colnames(DATA@meta.data) <- sub("SNN_res.","louvain_",colnames(DATA@meta.data))
  
  # for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
  #   temp2 <- DimPlot(DATA,dims = 1:2,reduction = i,group.by = sort(colnames(DATA@meta.data)[grep("louvain",colnames(DATA@meta.data))]), pt.size = .3,ncol = 8,label = T) +
  #     ggplot2::theme(legend.position = "none") 
  #   ggplot2::ggsave(temp2,filename = paste0("clustering_louvain_",i,".png"), path = paste0(opt$output_path,"/clustering"), dpi = 300,units = "mm",width = 140*8,height = 100*ceiling(length(grep("louvain",colnames(DATA@meta.data)))/8),limitsize = FALSE )
  # }
  
  for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
    s <- colnames(DATA@meta.data)[grep("louvain_",colnames(DATA@meta.data))]
    plot_list <- lapply(s,function(j){ DimPlot(DATA,dims = 1:2,reduction = i,group.by = j, pt.size = .3,ncol = 5 ,label = T) +
        ggplot2::theme(legend.position = "none") + ggplot2::ggtitle(label =j) })
    p <- cowplot::plot_grid(plotlist = plot_list, ncol=5)
    ggplot2::ggsave(p,filename = paste0("clustering_louvain_",i,".png"), path = paste0(opt$output_path,"/clustering"), dpi = 300,
                    units = "mm",width = 150*5,height = 140*ceiling(length(plot_list)/5),limitsize = FALSE )
  }
  
  }
rm(temp2); invisible(gc())
#---------



################################
### Clustering using Leiden ###
################################
if( 'leiden' %in% casefold(unlist(strsplit(opt$cluster_method,split = ","))) ){
  cat("\n### Clustering with leiden ###\n")
  if(!dir.exists(paste0(opt$output_path,"/clustering"))){dir.create(paste0(opt$output_path,"/clustering"))}
  #remove old clustering with this name
  DATA@meta.data <- DATA@meta.data[ , ! grepl( "leiden_", colnames(DATA@meta.data) ) ]
  
  DATA <- FindClusters(object = DATA, reduction.type = "pca", dims.use = 1:top_PCs, resolution = seq(.05,2,by=.05), verbose = T,graph.name = "SNN",algorithm = 4)
  colnames(DATA@meta.data) <- sub("SNN_res.","leiden_",colnames(DATA@meta.data))
  
  # for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
  #   temp2 <- DimPlot(DATA,dims = 1:2,reduction = i,group.by = sort(colnames(DATA@meta.data)[grep("leiden",colnames(DATA@meta.data))]), pt.size = .3,ncol = 8,label = T) +
  #     ggplot2::theme(legend.position = "none")
  #   ggplot2::ggsave(temp2,filename = paste0("clustering_leiden_",i,".png"), path = paste0(opt$output_path,"/clustering"), dpi = 300,units = "mm",width = 140*8,height = 100*ceiling(length(grep("leiden",colnames(DATA@meta.data)))/8),limitsize = FALSE )
  # }
  for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
    s <- colnames(DATA@meta.data)[grep("leiden_",colnames(DATA@meta.data))]
    plot_list <- lapply(s,function(j){ DimPlot(DATA,dims = 1:2,reduction = i,group.by = j, pt.size = .3,ncol = 5 ,label = T) +
        ggplot2::theme(legend.position = "none") + ggplot2::ggtitle(label = j) })
    p <- cowplot::plot_grid(plotlist = plot_list, ncol=5)
    ggplot2::ggsave(p,filename = paste0("clustering_leiden_",i,".png"), path = paste0(opt$output_path,"/clustering"), dpi = 300,
                    units = "mm",width = 150*5,height = 140*ceiling(length(plot_list)/5),limitsize = FALSE )
  }
  
  }
rm(temp2); invisible(gc())
#---------




######################################
### Hierachical clustering on UMAP ###
######################################
if( 'hc' %in% casefold(unlist(strsplit(opt$cluster_method,split = ","))) ){
  cat("\n### Clustering with HC (Hierachical CLustering on UMAP-10dims) ###\n")
  if(!dir.exists(paste0(opt$output_path,"/clustering"))){dir.create(paste0(opt$output_path,"/clustering"))}
  #remove old clustering with this name
  DATA@meta.data <- DATA@meta.data[ , ! grepl( "HC_", colnames(DATA@meta.data) ) ]
  
  h <- hclust(dist(DATA@reductions$umap10@cell.embeddings,method = "euclidean") ,method = "ward.D2")
  #h <- hclust(as.dist(DATA@graphs$SNN) ,method = "ward.D2")
  
  ideal <-round(sqrt(ncol(DATA)))
  for(k in sort(unique(round(ideal / seq(1,ideal,by = .3) ))) ){
    cat(k,"\t")  ;  cl <- cutree(h,k = k)
    DATA <- AddMetaData(DATA,metadata = setNames(cl,colnames(DATA)), paste0("HC_",k)) }

  for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
    s <- colnames(DATA@meta.data)[grep("HC_",colnames(DATA@meta.data))]
    plot_list <- lapply(s,function(j){ DimPlot(DATA,dims = 1:2,reduction = i,group.by = j, pt.size = .3,ncol = 5 ,label = T) +
        ggplot2::theme(legend.position = "none") + ggplot2::ggtitle(label = j) })
    p <- cowplot::plot_grid(plotlist = plot_list, ncol=5)
    ggplot2::ggsave(p,filename = paste0("clustering_HC_",i,".png"), path = paste0(opt$output_path,"/clustering"), dpi = 300,units = "mm",
                    width = 150*5,height = 140*ceiling(length(plot_list)/5),limitsize = FALSE )
  }
}
rm(temp2,h); invisible(gc())
#---------



####################################
### k-means partitioning on UMAP ###
####################################
if( 'kmeans' %in% casefold(unlist(strsplit(opt$cluster_method,split = ","))) ){
  if(!dir.exists(paste0(opt$output_path,"/clustering"))){dir.create(paste0(opt$output_path,"/clustering"))}
  #remove old clustering with this name
  DATA@meta.data <- DATA@meta.data[ , ! grepl( "kmeans_", colnames(DATA@meta.data) ) ]
  
  mincells <- 15
  ideal <- round(ncol(DATA) / mincells)
  clcl<- kmeans(DATA@reductions$umap10@cell.embeddings,centers = ncol(DATA)/10,iter.max = 50)
  DATA <- AddMetaData(DATA,metadata = setNames(clcl$cluster,colnames(DATA)), paste0("kmeans_",k))

  temp <- rowsum(DATA@reductions$umap10@cell.embeddings,clcl$cluster) / as.vector(table(clcl$cluster))
  cors <- cor(t(temp))

  cat("\nMerging clusters ...\n")
  merge_par <- seq(.8,.99,.02)

  for(j in merge_par){
    if( j > min(cors) ){
      tcors <- (cors > j)*1
      cell_clust <- clcl$cluster
      clust <- rownames(tcors)

      for( i in clust){
        sel <- rownames(tcors)[ tcors[i,] > 0 ]
        cell_clust[cell_clust %in% sel] <- sel[1]
      }
      DATA <- AddMetaData(object = DATA, metadata = cell_clust, col.name = paste0("kmeans_merged_",j))

      temp <- UMAPPlot(object = DATA, group.by=paste0("kmeans_merged_",j), pt.size = .5, plot.title= paste0("Clustering (kmeans_merged_",j,")"))+ggplot2::theme(legend.position = "none")
      ggsave(temp,filename = paste0("/clustering/UMAP_kmeans_merged_",j,".png"), path = opt$output_path, dpi = 300,units = "mm",width = 170,height = 150 )
    }
  }


  for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
    s <- colnames(DATA@meta.data)[grep("kmeans_",colnames(DATA@meta.data))]
    plot_list <- lapply(s,function(j){ DimPlot(DATA,dims = 1:2,reduction = i,group.by = j, pt.size = .3,ncol = 5 ,label = T) +
        ggplot2::theme(legend.position = "none") + ggplot2::ggtitle(label = j) })
    p <- cowplot::plot_grid(plotlist = plot_list, ncol=5)
    ggplot2::ggsave(p,filename = paste0("clustering_kmeans_",i,".png"), path = paste0(opt$output_path,"/clustering"), dpi = 300,
                    units = "mm",width = 150*5,height = 140*ceiling(length(plot_list)/5),limitsize = FALSE )
  }
}
#---------



#######################
### HBDSCAN on UMAP ###
#######################
if( 'hdbscan' %in% casefold(unlist(strsplit(opt$cluster_method,split = ","))) ){
  cat("\n### Clustering with HDBSCAN on UMAP-10dims ###\n")
  if(!dir.exists(paste0(opt$output_path,"/clustering"))){dir.create(paste0(opt$output_path,"/clustering"))}
  #remove old clustering with this name
  DATA@meta.data <- DATA@meta.data[ , ! grepl( "hdbscan_", colnames(DATA@meta.data) ) ]
  
  for(k in seq(5,100,by=2)){
    cat(k,"\t")
    clusters <- hdbscan(DATA@reductions$umap@cell.embeddings, minPts = k)
    names(clusters$cluster) <- rownames(DATA@meta.data)
    DATA <- AddMetaData(object = DATA, metadata = clusters$cluster, col.name = paste0("hdbscan_",k))
  }
  # for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
  #   temp2 <- DimPlot(DATA,dims = 1:2,reduction = i,group.by = colnames(DATA@meta.data)[grep("hdbscan_",colnames(DATA@meta.data))], pt.size = .3,ncol = 5)
  #   ggplot2::ggsave(temp2,filename = paste0("clustering_hdbscan_",i,".png"), path = paste0(opt$output_path,"/clustering"), dpi = 300,
  #                   units = "mm",width = 140*5,height = 100*ceiling(length(grep("hdbscan_",colnames(DATA@meta.data)))/5),limitsize = FALSE )
  # }
  
  for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
    s <- colnames(DATA@meta.data)[grep("hdbscan_",colnames(DATA@meta.data))]
    plot_list <- lapply(s,function(j){ DimPlot(DATA,dims = 1:2,reduction = i,group.by = j, pt.size = .3,ncol = 5 ,label = T) +
        ggplot2::theme(legend.position = "none") + ggplot2::ggtitle(label = j) })
    p <- cowplot::plot_grid(plotlist = plot_list, ncol=5)
    ggplot2::ggsave(p,filename = paste0("clustering_hdbscan_",i,".png"), path = paste0(opt$output_path,"/clustering"), dpi = 300,
                    units = "mm",width = 150*5,height = 140*ceiling(length(plot_list)/5),limitsize = FALSE )
  }
  
}
rm(temp2); invisible(gc())
#---------



###################################
### SAVING RAW Seurat.v3 OBJECT ###
###################################
cat("\n### Saving the Seurat object ###\n")
saveRDS(DATA, file = paste0(opt$output_path,"/seurat_object.rds") )
write.csv2(DATA@meta.data,paste0(opt$output_path,"/Metadata_with_clustering.csv"),row.names = T)
#---------



#############################
### SYSTEM & SESSION INFO ###
#############################
#---------
print_session_info()
#---------
