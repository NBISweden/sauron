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
  make_option(c("-p", "--PCs_use"),               type = "character",   metavar="character",   default='var,1', help="Method and threshold level for selection of significant principal components. The method should be separated from the threshold via a comma. 'top,5' will use the top 5 PCs, which is the default. 'var,1' will use all PCs with variance above 1%."),
  make_option(c("-v", "--var_genes"),             type = "character",   metavar="character",   default='scran',  help="Whether use 'Seurat' or the 'Scran' method for variable genes identification. An additional value can be placed after a comma to define the level of dispersion wanted for variable gene selection. 'Seurat,2' will use the threshold 2 for gene dispersions. Defult is 'scran'. For Scran, the user should inpup the level of biological variance 'Scran,0.2'. An additional blocking parameter (a column from the metadata) can ba supplied to 'Scran' method block variation comming from uninteresting factors, which can be parsed as 'Scran,0.2,Batch'."),
  make_option(c("-s", "--cluster_use"),           type = "character",   metavar="character",   default='all',    help="The cluster of cells to be used for analysis. Should be defined as the clustering name followed by the cluster names to be used, comma-separated. E.g.: 'SNN_0.2,1,2,3,5,6'."),
  make_option(c("-m", "--cluster_method"),        type = "character",   metavar="character",   default='snn,hc', help="The clustering method and cluster to select for analysis. Current methods are 'hc','snn','dbscan','hdbscan','flowpeaks'. If no input is suplied, all methods will be run."),
  make_option(c("-d", "--dim_reduct_use"),        type = "character",   metavar="character",   default='umap',  help="Which dimensionality reduction method to be run on top of PCA: UMAP (default) or tSNE. If both, then specify them comma-separated'UMAP,tSNE'."),
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
pkgs <- c("Seurat","rafalib","scran","biomaRt","scater","dplyr","RColorBrewer","dbscan","flowPeaks","scales","igraph","sva","parallel")
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
    cells_use <- rownames(DATA@meta.data)[factor(DATA@meta.data[,clustering_use]) %in% clusters_to_select]}   #Filter out cells with no assigned clusters
  DATA <- SubsetData(DATA,assay = DefaultAssay(DATA),cells = cells_use)
} else {
  cat("\nThe name of the cluster or the cluster name were not found in your data.\n All cells will be used ...\n")
  cells_use <- rownames(DATA@meta.data)}
# sel <- rowSums(as.matrix(DATA@assays$RNA@counts) >= 1) >= 1
# DATA <- CreateSeuratObject(as.matrix(DATA@assays$RNA@counts[sel,cells_use]), meta.data = DATA@meta.data[cells_use,])
#---------



###########################
### FIND VARIABLE GENES ###
###########################
if(DefaultAssay(DATA) == "RNA"){
  output_path <- paste0(opt$output_path,"/variable_genes")
  DATA <- compute_hvgs(DATA,VAR_choice,output_path)
} else { DATA@assays[[DefaultAssay(DATA)]]@var.features <- rownames(DATA@assays[[DefaultAssay(DATA)]]@data)}
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
}
#---------



###################
### Running PCA ###
###################
cat("\n### Running PCA ###\n")
if(!dir.exists(paste0(opt$output_path,"/pca_plots"))){dir.create(paste0(opt$output_path,"/pca_plots"),recursive = T)}
DATA <- RunPCA(DATA, do.print = F, pcs.compute = 50, assay = DefaultAssay(DATA))
write.csv2(DATA@reductions$pca@cell.embeddings, paste0(opt$output_path,"/pca_plots/PCA_coordinates.csv"))
write.csv2(DATA@reductions$pca@feature.loadings, paste0(opt$output_path,"/pca_plots/PCA_feature_loadings.csv"))

ggsave(PCHeatmap(DATA,ncol=5,dims=1:10),filename = paste0("PCA_heatmap.png"), path = paste0(opt$output_path,"/pca_plots"), dpi = 300,units = "mm",width = 150*5,height = 150*4.5)

var_expl <- (DATA@reductions$pca@stdev^2)/sum(DATA@reductions$pca@stdev^2)
PC_choice <- as.character(unlist(strsplit(opt$PCs_use,",")))
if(casefold(PC_choice[1]) == "var"){  top_PCs <- sum( var_expl > as.numeric(PC_choice[2])/100 )
} else if(casefold(PC_choice[1]) == "top"){top_PCs <- as.numeric(PC_choice[2])}

png(filename = paste0(opt$output_path,"/pca_plots/Varianc_explained_PC.png"),width = 1500,height =1200,res = 200)
plot( var_expl,yaxs="i",bg="grey",pch=21,type="l",ylab="% Variance",xlab="PCs",main="all cells",las=1)
points( var_expl,bg=c(rep("orange",top_PCs),rep("grey",50-top_PCs)),pch=21)
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
    DATA@reductions[["tsne"]] <- CreateDimReducObject(embeddings = as.matrix(read.csv2(paste0(opt$output_path,"/tsne_plots/tSNE_coordinates.csv"),row.names = 1)),key = "tSNE_",assay = DefaultAssay(DATA))
  } else { cat("\nPre-computed tSNE NOT found. Computing tSNE ...\n")
    if(DefaultAssay(DATA) == "RNA"){n <- 1:top_PCs } else { n <- 1:50 }
    ttt <- Sys.time()
    DATA <- RunTSNE(object = DATA, perplexity=50, max_iter=2000,theta=0.1,eta=2000,exaggeration_factor=12,dims.use = n,verbose = T,num_threads=0)
    cat("multicore tSNE ran in ",difftime(Sys.time(), ttt, units='mins'))
    write.csv2(DATA@reductions$tsne@cell.embeddings, paste0(opt$output_path,"/tsne_plots/tSNE_coordinates.csv"))}
}
#---------



####################
### Running UMAP ###
####################
if( "umap" %in% casefold(unlist(strsplit(opt$dim_reduct_use,",")))){
  cat("\n### Running UMAP ###\n")
  if(!dir.exists(paste0(opt$output_path,"/_plots"))){dir.create(paste0(opt$output_path,"/umap_plots"),recursive = T)}

  if(file.exists(paste0(opt$output_path,"/umap_plots/UMAP_coordinates.csv"))&file.exists(paste0(opt$output_path,"/umap_plots/UMAP10_coordinates.csv"))){
    cat("\nPre-computed UMAP found and will be used:\n",paste0(opt$output_path,"/umap_plots/UMAP_coordinates.csv"),"\n")
    DATA@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(read.csv2(paste0(opt$output_path,"/umap_plots/UMAP_coordinates.csv"),row.names = 1)),key = "UMAP_",assay = DefaultAssay(DATA))
    DATA@reductions[["umap10"]] <- CreateDimReducObject(embeddings = as.matrix(read.csv2(paste0(opt$output_path,"/umap_plots/UMAP10_coordinates.csv"),row.names = 1)),key = "UMAP_",assay = DefaultAssay(DATA))

  } else { cat("\nPre-computed UMAP NOT found. Computing UMAP ...\n")
    if(DefaultAssay(DATA) == "RNA"){n <- 1:top_PCs } else { n <- 1:50 }
    ttt <- Sys.time()
    DATA <- RunUMAP(object = DATA, dims = n,n.components = 2, n.neighbors = 50,min.dist = 0.0001, verbose = T,num_threads=0)
    cat("UMAP_2dimensions ran in ",difftime(Sys.time(), ttt, units='mins'))
    invisible(gc())
    ttt <- Sys.time()
    
    DATA <- RunUMAP(object = DATA, dims = n,n.components = 10, n.neighbors = 10, verbose = T,num_threads=0,reduction.name = "umap10",reduction.key = "umap10_")
    cat("UMAP_10dimensions ran in ",difftime(Sys.time(), ttt, units='mins'))
    invisible(gc())
    write.csv2(DATA@reductions$umap@cell.embeddings, paste0(opt$output_path,"/umap_plots/UMAP_coordinates.csv"))
    write.csv2(DATA@reductions$umap10@cell.embeddings, paste0(opt$output_path,"/umap_plots/UMAP10_coordinates.csv"))}
  invisible(gc())
}
#---------



###################
### Running SNN ###
###################
cat("\n### Running SNN ###\n")
DATA <- FindNeighbors(DATA,assay = DefaultAssay(DATA),graph.name="SNN")
g <- graph_from_adjacency_matrix(as.matrix(DATA@graphs$SNN),weighted = T,diag = F,mode = "undirected")
saveRDS(DATA@graphs$SNN, file = paste0(opt$output_path,"/SNN_Graph.rds") )
for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
  png(filename = paste0(opt$output_path,"/",i,"_plots","/",i,"_plot_with_SNN_overlay.png"),width = 200,height =205,res = 600,units = "mm")
  plot.igraph(x = g, layout = DATA@reductions[[i]]@cell.embeddings[,1:2], edge.width = E(graph = g)$weight/4, vertex.label = NA,
              edge.color = colorRampPalette(c("grey90","black"))(50)[round(E(graph = g)$weight/4*49+1)],
              vertex.size = 1,vertex.frame.color=hue_pal(l=50, c=80)(length(levels(factor(DATA$orig.ident))))[factor(DATA$orig.ident)],
              vertex.color = hue_pal()(length(unique(DATA$orig.ident)))[factor(DATA$orig.ident)])
  invisible(dev.off())
}
rm(g); invisible(gc())
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
rm(temp,temp2); invisible(gc())
#---------



################################
### Clustering using Louvain ###
################################
if( 'snn' %in% casefold(unlist(strsplit(opt$cluster_method,split = ","))) ){
  cat("\n### Clustering with SNN ###\n")
  if(!dir.exists(paste0(opt$output_path,"/clustering"))){dir.create(paste0(opt$output_path,"/clustering"))}
  for(k in seq(.05,2,by=.05)){
    cat(k,"\t")
    DATA <- FindClusters(object = DATA, reduction.type = "pca", dims.use = 1:top_PCs, resolution = k, verbose = F,graph.name = "SNN")
  }
  for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
    temp2 <- DimPlot(DATA,dims = 1:2,reduction = i,group.by = sort(colnames(DATA@meta.data)[grep("SNN",colnames(DATA@meta.data))]), pt.size = .3,ncol = 8)
    ggplot2::ggsave(temp2,filename = paste0("clustering_SNN_",i,".png"), path = paste0(opt$output_path,"/clustering"), dpi = 300,units = "mm",width = 140*8,height = 100*ceiling(length(grep("SNN",colnames(DATA@meta.data)))/8),limitsize = FALSE )
  }}
rm(temp2); invisible(gc())
#---------



######################################
### Hierachical clustering on UMAP ###
######################################
if( 'hc' %in% casefold(unlist(strsplit(opt$cluster_method,split = ","))) ){
  cat("\n### Clustering with HC (Hierachical CLustering on UMAP-10dims) ###\n")
  if(!dir.exists(paste0(opt$output_path,"/clustering"))){dir.create(paste0(opt$output_path,"/clustering"))}
  h <- hclust(dist(DATA@reductions$umap10@cell.embeddings,method = "euclidean") ,method = "ward.D2")

  ideal <-round(sqrt(ncol(DATA)))
  for(k in sort(unique(round(ideal / seq(1,ideal,by = .3) ))) ){
    cat(k,"\t")  ;  cl <- cutree(h,k = k)
    DATA <- AddMetaData(DATA,metadata = setNames(cl,colnames(DATA)), paste0("HC_",k)) }

  for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
    s <- colnames(DATA@meta.data)[grep("HC_",colnames(DATA@meta.data))]
    plot_list <- lapply(s,function(j){ DimPlot(DATA,dims = 1:2,reduction = i,group.by = j, pt.size = .3,ncol = 8 ,label = T) +
        ggplot2::theme(legend.position = "none") + ggplot2::ggtitle(label = paste0("Hierc.Clust (",j,")")) })
    p <- cowplot::plot_grid(plotlist = plot_list, ncol=8)
    ggplot2::ggsave(p,filename = paste0("clustering_HC_",i,".png"), path = paste0(opt$output_path,"/clustering"), dpi = 300,units = "mm",width = 150*8,height = 140*ceiling(length(plot_list)/8),limitsize = FALSE )
  }
}
rm(temp2,h); invisible(gc())
#---------



####################################
### k-means partitioning on UMAP ###
####################################
if( 'kmeans' %in% casefold(unlist(strsplit(opt$cluster_method,split = ","))) ){
  if(!dir.exists(paste0(opt$output_path,"/clustering"))){dir.create(paste0(opt$output_path,"/clustering"))}

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
    plot_list <- lapply(s,function(j){ DimPlot(DATA,dims = 1:2,reduction = i,group.by = j, pt.size = .3,ncol = 8 ,label = T) +
        ggplot2::theme(legend.position = "none") + ggplot2::ggtitle(label = paste0("Kmeans.Clust (",j,")")) })
    p <- cowplot::plot_grid(plotlist = plot_list, ncol=8)
    ggplot2::ggsave(p,filename = paste0("clustering_kmeans_",i,".png"), path = paste0(opt$output_path,"/clustering"), dpi = 300,units = "mm",width = 150*8,height = 140*ceiling(length(plot_list)/8),limitsize = FALSE )
  }
}
#---------



#######################
### HBDSCAN on UMAP ###
#######################
if( 'hdbscan' %in% casefold(unlist(strsplit(opt$cluster_method,split = ","))) ){
  cat("\n### Clustering with HDBSCAN on UMAP-10dims ###\n")
  if(!dir.exists(paste0(opt$output_path,"/clustering"))){dir.create(paste0(opt$output_path,"/clustering"))}
  for(k in seq(5,100,by=2)){
    cat(k,"\t")
    clusters <- hdbscan(DATA@reductions$umap@cell.embeddings, minPts = k)
    names(clusters$cluster) <- rownames(DATA@meta.data)
    DATA <- AddMetaData(object = DATA, metadata = clusters$cluster, col.name = paste0("hdbscan_",k))
  }
  for(i in c("pca",casefold(unlist(strsplit(opt$dim_reduct_use,","))))){
    temp2 <- DimPlot(DATA,dims = 1:2,reduction = i,group.by = colnames(DATA@meta.data)[grep("hdbscan_",colnames(DATA@meta.data))], pt.size = .3,ncol = 8)
    ggplot2::ggsave(temp2,filename = paste0("clustering_hdbscan_",i,".png"), path = paste0(opt$output_path,"/clustering"), dpi = 300,units = "mm",width = 140*8,height = 100*ceiling(length(grep("hdbscan_",colnames(DATA@meta.data)))/8),limitsize = FALSE )
  }
}
rm(temp2); invisible(gc())
#---------



###################################
### SAVING RAW Seurat.v3 OBJECT ###
###################################
cat("\n### Saving the Seurat object ###\n")
saveRDS(DATA, file = paste0(opt$output_path,"/Seurat_object.rds") )
write.csv2(DATA@meta.data,paste0(opt$output_path,"/Metadata_with_clustering.csv"),row.names = T)
#---------



#############################
### SYSTEM & SESSION INFO ###
#############################
#---------
print_session_info()
#---------
