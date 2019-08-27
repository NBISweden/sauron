#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
cat("\nRunning TRAJECTORY ANALYSIS with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object FILE."),
  make_option(c("-c", "--clustering_use"),        type = "character",   metavar="character",   default='none',  help="The clustering column to be used. Should be chosen from one of the method in script 02."),
  make_option(c("-m", "--metadata_use"),          type = "character",   metavar="character",   default='none',  help="Column names of the metadata table to plot in the trajectory map."),
  make_option(c("-d", "--reduction_use"),         type = "character",   metavar="character",   default='umap10',  help="Dimensionality reduction method to base the trajectories on. It could be a pre-computed slot within your Seurat Object (in case: pca, umap, umap10, tsne) or it will compute additional others (ica, dm). Trajectories are defined on top of UMAP ambedding with 10 dimensions by default."),
  make_option(c("-r", "--reduction_visualize"),   type = "character",   metavar="character",   default='umap',  help="Dimensionality reduction method to visualize trajectories. It could be a pre-computed slot within your Seurat Object (in case: pca, umap, umap10, tsne) or it will compute additional others (ica, dm). Trajectories are defined on top of UMAP ambedding with 10 dimensions by default."),
  make_option(c("-e", "--method_use"),            type = "character",   metavar="character",   default='DDRTree',  help="Method used for trajectory estimation, MST or DDRTree. DDRTree is default."),
  make_option(c("-t", "--no_traj_components"),    type = "character",   metavar="character",   default='3',  help="Number trajectory components to use"),
  make_option(c("-p", "--no_of_paths"),           type = "character",   metavar="character",   default='none',  help="Number specifying how many paths to allow."),
  make_option(c("-s", "--start_cluster"),         type = "character",   metavar="character",   default='none',  help="Cluster from which the trajectories will start from."),
  make_option(c("-a", "--assay"),                 type = "character",   metavar="character",   default='RNA',  help="Assay to use for trajectory differential expression. Default is 'RNA'"),
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
pkgs <- c("Seurat","scales","monocle","plot3D","scales","parallel")
inst_packages(pkgs)
#---------



#############################
### LOAD Seurat.v3 OBJECT ###
#############################
cat("\n### LOADING Seurat.v3 OBJECT ###\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------



#########################
### DEFINE PARAMETERS ###
#########################
reduction_use <- unlist(strsplit( opt$reduction_use , split = ","))
red_use <- unlist(strsplit( opt$reduction_use , split = ","))
red_vis <- unlist(strsplit( opt$reduction_visualize , split = ","))
#---------



########################
### BUILD SCE OBJECT ###
########################
cat("\nConstructing monocle object ...\n")
#For the whole dataset as monocle object
temp <- newCellDataSet(cellData = as.matrix(DATA@assays[[opt$assay]]@counts[rownames(DATA@assays[[opt$assay]]@counts),]), 
        phenoData = new("AnnotatedDataFrame", data = cbind(DATA@meta.data ) ) ,
        featureData = new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(DATA@assays[[opt$assay]]@counts),row.names = rownames(DATA@assays[[opt$assay]]@counts))) )
temp <- estimateSizeFactors(temp)

#For the whole data
if( casefold(opt$method_use) == "ddrtree"){
  temp2 <- newCellDataSet(cellData = t(DATA@reductions[[opt$reduction_use]]@cell.embeddings),
                          phenoData = new("AnnotatedDataFrame", data = cbind(DATA@meta.data ) ) )
  temp2 <- estimateSizeFactors(temp2)
  tempUMAP2 <- reduceDimension(temp2, reduction_method = "DDRTree",norm_method = "none",max_components = as.numeric(opt$no_traj_components))
} else if(casefold(opt$method_use) == "ica"){
  tempUMAP2 <- reduceDimension( temp, reduction_method = "ICA",norm_method = "none",max_components = as.numeric(opt$no_traj_components))
  tempUMAP2@reducedDimS <- t(opt$reduction_use)
}
#---------



###################
### ORDER CELLS ###
###################
cat("\nOrdering cells ...\n")
tempUMAP2 <- orderCells(tempUMAP2)

if( ( opt$no_of_paths != "none" ) & ( opt$start_cluster != "none" ) ){
  tempUMAP2 <- orderCells(tempUMAP2 , root_state = as.numeric(opt$start_cluster),num_paths = as.numeric(opt$no_of_paths))
} else if( opt$no_of_paths != "none" ){
  tempUMAP2 <- orderCells(tempUMAP2 , num_paths = as.numeric(opt$no_of_paths))
} else if( opt$start_cluster != "none" ){
  tempUMAP2 <- orderCells(tempUMAP2 , root_state = as.numeric(opt$start_cluster))
}


cat("\nReassigning gene expression to object \n")
if(file.exists(paste0(opt$output_path,"/Monocle_trajectory_object.rds"))){
  tempUMAP3<- readRDS(paste0(opt$output_path,"/Monocle_trajectory_object.rds"))
} else {
  tempUMAP3 <- tempUMAP2
  tempUMAP3@assayData$exprs <- temp@assayData$exprs
  tempUMAP3@featureData <- temp@featureData
  tempUMAP3 <- estimateSizeFactors(tempUMAP3)

  saveRDS(tempUMAP3,paste0(opt$output_path,"/Monocle_trajectory_object.rds"))
}
#---------



########################
### PLOTTING RESULTS ###
########################
cat("\nPlotting results ...\n")

for(i in c(opt$clustering_use,unlist(strsplit(opt$metadata_use,",")),"State","Pseudotime") ){
  message(i)
  pdf(paste0(opt$output_path,"/Trajectories_by_",i,".pdf"),width = 8,height = 8,useDingbats = F)
  print(plot_cell_trajectory(tempUMAP3, color_by = i,y=2,x=1,cell_size = .5))
  
  print(plot_complex_cell_trajectory(tempUMAP3, color_by = i, show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3))
 
  x <- factor(tempUMAP3@phenoData@data[,i])
  scatter3D(tempUMAP3@reducedDimS[1,],tempUMAP3@reducedDimS[2,],tempUMAP3@reducedDimS[3,], bty = "g", theta = 40,
            colvar = as.integer(factor(tempUMAP3@phenoData@data[,i])), col =  hue_pal()(length(levels(x))),
            pch=16,cex=.5,xlab = "Dim1", ylab ="Dim2", zlab = "Dim3",colkey = list(at = seq(1.5,length(levels(x))+.5,length.out = length(levels(x))), 
            side = 1, labels = levels(x) ) )
  scatter3D(tempUMAP3@reducedDimK[1,],tempUMAP3@reducedDimK[2,],tempUMAP3@reducedDimK[3,],add=T,col="black",pch=16,cex=.5)
  
  for(j in names(DATA@reductions) ){
    message("   ",j)
    plot(DATA@reductions[[j]]@cell.embeddings[,1:2],col=hue_pal()(length(levels(x)))[x],cex=.7,pch=16,main=i,xlab = paste0(j,"_1"),ylab = paste0(j,"_2"),las=1)
  }

dev.off()
}
#---------




#########################################
### COMPUTING DGE PER BRANCHING POINT ###
#########################################
cat("\nComputing DGE per Branch ...\n")
n_branches <- length(tempUMAP3@auxOrderingData$DDRTree$branch_points)

for(i in 1:n_branches){
  message('processing branch ',i)
  if(file.exists(paste0(opt$output_path,"/diff_test_branch_",i,".csv"))){
    BEAM_res <- read.csv2(paste0(opt$output_path,"/diff_test_branch_",i,".csv"),row.names = 1)
  } else {
    BEAM_res <- BEAM(tempUMAP3, branch_point = i, cores = detectCores()-1)
    BEAM_res <- BEAM_res[order(BEAM_res$qval),]
    write.csv2(BEAM_res,paste0(opt$output_path,"/diff_test_branch_",i,".csv"),row.names = T)
  }
  n <- min(100, length(row.names(subset(BEAM_res, qval < 1e-6))))
  pdf(paste0(opt$output_path,"/Trajectory_branch",i,".pdf"),width = 5,height = 6,useDingbats = F)
  try(plot_genes_branched_heatmap(tempUMAP3[rownames(BEAM_res)[1:n],],branch_point = i, norm_method= "log",
                                  num_clusters = 4, use_gene_short_name = T,  show_rownames = T))
  dev.off()
}

#---------



###################################
### SAVING RAW Seurat.v3 OBJECT ###
###################################
cat("\n### Saving the updated back to the Seurat object used as input###\n")
DATA$Pseudotime <- tempUMAP3$Pseudotime
DATA$State <- tempUMAP3$State
DATA@reductions[[opt$method_use]] <- CreateDimReducObject(embeddings = as.matrix(tempUMAP3@reducedDimS), key = opt$method_use, assay = opt$assay)

saveRDS(DATA, file = opt$Seurat_object_path )
write.csv2(DATA@meta.data,paste0(opt$output_path,"/Metadata_with_trajectory.csv"),row.names = T)
#---------



#############################
### SYSTEM & SESSION INFO ###
#############################
#---------
print_session_info()
#---------
