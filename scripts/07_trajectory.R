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
  make_option(c("-m", "--metadata_use"),          type = "character",   metavar="character",   default='none',  help="Column names of the metadata table to plot in the trajectory map."),
  make_option(c("-e", "--method_use"),            type = "character",   metavar="character",   default='none',  help="Method used for trajectory estimation, MST or DDRTree. Trajectories are defined on top of UMAP ambedding with 10 dimensions."),
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
pkgs <- c("rafalib","dplyr","RColorBrewer","scales","igraph","monocle","Seurat","pheatmap","plot3D")
inst_packages(pkgs)
#---------




### LOAD Seurat OBJECT 
#---------
cat("\nLoading/ data and metadata ...\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------





### Run UMAP with more dimensions
#---------
if( length(DATA@dr$umap) != 1){
  cat("\nUMAP not found in the dataset ... COmputing UMAP with 2 and 10 dimensions\n")
  DATA <- RunUMAP(DATA,verbose=T,max.dim = 2,min_dist = 1e-5,n_neighbors = 50,dims.use = 1:20)
  DATA <- RunUMAP(DATA,verbose=T,max.dim = 10,min_dist = 1e-5,n_neighbors = 50,dims.use = 1:20,reduction.name="umap10")
}

umap2 <- DATA@dr$umap@cell.embeddings
umap10 <- DATA@dr$umap10@cell.embeddings
#---------





### Constructing Monocle object
#---------
cat("\nConstructing monocle object ...\n")

#For the whole dataset as monocle object
temp <- newCellDataSet(cellData = as.matrix(DATA@data[rownames(DATA@data),]), 
        phenoData = new("AnnotatedDataFrame", data = cbind(DATA@meta.data ) ) ,
        featureData = new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(DATA@data),row.names = rownames(DATA@data))) )
temp <- estimateSizeFactors(temp)

#For the whole data
if( casefold(opt$method_use) == "ddrtree"){
  temp2 <- newCellDataSet(cellData = t(umap10), phenoData = new("AnnotatedDataFrame", data = cbind(DATA@meta.data ) ) )
  temp2 <- estimateSizeFactors(temp2)
  tempUMAP2 <- reduceDimension(temp2,reduction_method = "DDRTree",norm_method = "none",max_components = 3)
} else if(casefold(opt$method_use) == "ica"){
  tempUMAP2 <- reduceDimension(temp,reduction_method = "ICA",norm_method = "none",max_components = 3)
  tempUMAP2@reducedDimS <- t(umap10)
}
#---------




### Constructing Monocle object
#---------
cat("\nOrdering cells ...\n")

tempUMAP2 <- orderCells(tempUMAP2)

pdf(paste0(output_path,"/Trajectory_",opt$method_use,"_default.pdf"),width = 8,height = 8,useDingbats = F)
for(i in c(opt$metadata_use,"State","Pseudotime") ){
  print(plot_cell_trajectory(tempUMAP2, color_by = i,y=2,x=1,cell_size = .5))
}

scatter3D(tempUMAP2@reducedDimS[1,],tempUMAP2@reducedDimS[2,],tempUMAP2@reducedDimS[3,], bty = "g", theta = 20,
          colvar = as.integer(factor(tempUMAP2@phenoData@data$res.0.4,levels = 0:6)), col =  hue_pal()(7),
          pch=16,cex=.5,xlab = "Dim1", ylab ="Dim2", zlab = "Dim3",colkey = list(at = seq(1.5,6.5,length.out = length(levels(factor(tempUMAP2@phenoData@data$res.0.4)))), 
          side = 1, labels = levels(factor(tempUMAP2@phenoData@data$res.0.4,levels = 0:6)) ) )
scatter3D(tempUMAP2@reducedDimK[1,],tempUMAP2@reducedDimK[2,],tempUMAP2@reducedDimK[3,],add=T,col="black",pch=16,cex=.5)

x <- factor(tempUMAP2@phenoData@data[,"State"])
scatter3D(tempUMAP2@reducedDimS[1,],tempUMAP2@reducedDimS[2,],tempUMAP2@reducedDimS[3,], bty = "g", theta = 40,
          colvar = as.integer(factor(tempUMAP2@phenoData@data$State)), col =  hue_pal()(length(levels(x))),
          pch=16,cex=.5,xlab = "Dim1", ylab ="Dim2", zlab = "Dim3",colkey = list(at = seq(1.5,length(levels(x))+.5,length.out = length(levels(x))), 
          side = 1, labels = levels(factor(tempUMAP2@phenoData@data$State)) ) )
scatter3D(tempUMAP2@reducedDimK[1,],tempUMAP2@reducedDimK[2,],tempUMAP2@reducedDimK[3,],add=T,col="black",pch=16,cex=.5)


plot(umap2[,],col=hue_pal()(length(levels(tempUMAP2$State)))[tempUMAP2$State],cex=.7,pch=16,main="State")
plot(umap2[,],col=hue_pal()(length(levels(factor(tempUMAP2$res.0.4))))[factor(tempUMAP2$res.0.4)],cex=.7,pch=16,main="res.0.4")
plot(umap2[,],col=hue_pal()(length(levels(factor(tempUMAP2$days_post_infection))))[factor(tempUMAP2$days_post_infection)],cex=.7,pch=16,main="days_post_infection")
dev.off()
#---------








### Constructing Monocle object
#---------
cat("\nOrdering cells ...\n")

tempUMAP2@

#---------




























