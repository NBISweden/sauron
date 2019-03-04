#!/usr/bin/env Rscript

### LOAD OPTPARSE
#---------
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')};
library(optparse)
#---------


### DEFINE PATH TO LOCAL FILES
#---------
cat("\nRunning CELL TYPE PREDICTION with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object"),
  make_option(c("-m", "--marker_lists"),          type = "character",   metavar="character",   default='none',  help="A folder containing .csv files with the list of markers for comparison"),
  make_option(c("-c", "--cluster_use"),           type = "character",   metavar="character",   default='all',   help="The clustering name to be used for the analysis"),
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
#if(!dir.exists(paste0(opt$output_path,"/cell_prediction"))){dir.create(paste0(opt$output_path,"/cell_prediction"),recursive = T)}
#---------



### LOAD LIBRARIES
#---------
cat("\nLoading/installing libraries ...\n")
initial.options <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(sub("--file=","",initial.options[grep("--file=",initial.options)]))
source( paste0(script_path,"/inst_packages.R") )
pkgs <- c("Seurat","rafalib","scran","biomaRt","scater","dplyr","RColorBrewer","dbscan","flowPeaks","scales","igraph","sva")
inst_packages(pkgs)
#---------




### LOAD Seurat OBJECT 
#---------
cat("\nLoading/ data and metadata ...\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------






### LOAD cell marker list
#---------
cat("\nLoading cell marker lists and creating cell identity matrix ...\n")
cat("\nThe following datasets were found ...\n")
marker_lists <- list.files( opt$marker_lists , pattern = ".csv")
print(marker_lists)

for(i in marker_lists){

cat("\nProcessing list ", i ," ...\n")
cellIDs <- read.csv2(paste0(opt$marker_lists,"/",i),header = F)
rownames(cellIDs) <- as.character(cellIDs[,1])
print(cellIDs)
cellIDs <- as.list(as.data.frame(t(cellIDs[,-1])))
cellIDs <- lapply(cellIDs, function(x) as.character(x[x!=""]) )




#Create cell identity matrix
cell_ident <- unique(unlist(cellIDs))
#cell_ident <- lapply(cellIDs,function(x) cell_ident %in% x)
cell_ident <- lapply(cellIDs,function(x) ifelse(cell_ident %in% x,1,ifelse(paste0("-",cell_ident) %in% x,-1,0)))
cell_ident <- as.data.frame(cell_ident)
rownames(cell_ident) <- unique(unlist(cellIDs))
#--------------------------------



#Filter data for the genes inthe cell identity matrix
#--------------------------------
data <- as.matrix(DATA@data)
sel <- rownames(cell_ident)[ rownames(cell_ident) %in% rownames(data) ]
data <- data[sel,]
cell_ident <- cell_ident[sel,]
gc()
#--------------------------------

#Using correlations as distance
#--------------------------------
cat("\nPredicted cell type by correlation method ...\n")
cors <- apply(data,2,function(x) cor(x , cell_ident) )
rownames(cors) <- colnames(cell_ident)
cors <- t(t(cors) / apply(cors,2,max))
cors[1:nrow(cors),1:20]

write.csv(cors,paste0(opt$output_path,"/Cell_pred_correlation_",i,".csv"),row.names = T)

pred <- unlist( apply(cors,2,function(x) colnames(cell_ident) [which.max(x)]) )
my_nas <- colnames(cors)[! colnames(cors) %in% names(pred)]
pred <- c(pred , setNames(rep(NA,length(my_nas)),my_nas))

print(table(pred))
DATA <- AddMetaData(DATA,metadata = pred,col.name = paste0("cell_pred_correlation_",i))

png(filename = paste0(opt$output_path,"/tSNE_cell_pred_correlation_",i,".png"),width = 700,height = 600,res = 150)
TSNEPlot(object = DATA,group.by=paste0("cell_pred_correlation_",i),pt.size = .3,plot.title= paste0("cell_pred_correlation_",i))
invisible(dev.off())
#--------------------------------


#Using normalized euclidean distance
#--------------------------------
cat("\nPredicted cell type by euclidean-distance method ...\n")
c <- apply(data,2,function(x) cell_ident - (x/max(x)) )
c1 <- lapply(c,function(x) colSums(x^2) )
c2 <- as.data.frame(c1)
c2[1:nrow(c2),1:10]

write.csv(c2,paste0(opt$output_path,"/cell_pred_euclidean_",i,".csv"),row.names = T)


pred2 <- unlist(apply(c2,2,function(x) colnames(cell_ident) [x==min(x)]))
my_nas <- colnames(cors)[! colnames(cors) %in% names(pred2)]
pred2 <- c(pred2 , setNames(rep(NA,length(my_nas)),my_nas))

print(table(pred2))
DATA <- AddMetaData(DATA,metadata = pred2,col.name = paste0("cell_pred_euclidean_",i))
png(filename = paste0(opt$output_path,"/tSNE_cell_pred_euclidean_",i,".png"),width = 700,height = 600,res = 150)
TSNEPlot(object = DATA,group.by=paste0("cell_pred_euclidean_",i),pt.size = .3,plot.title= paste0("cell_pred_euclidean_",i))
invisible(dev.off())
#--------------------------------
}





### Saving the Seurat object
#---------
saveRDS(DATA, file = opt$Seurat_object_path )
#---------





cat("\n!!! Script executed Sucessfully !!!\n")


### System and session information
#---------
cat("\n\n\n\n... SYSTEM INFORMATION ...\n")
INFORMATION <- Sys.info()
print(as.data.frame(INFORMATION))

cat("\n\n\n\n... SESSION INFORMATION ...\n")
sessionInfo()
#---------



