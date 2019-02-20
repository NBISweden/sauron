#!/usr/bin/env Rscript

### LOAD OPTPARSE
#---------
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')};
library(optparse)
#---------


### DEFINE PATH TO LOCAL FILES
#---------
cat("\nCREATING SEURAT OBJECT with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--input_path"),       type = "character",   metavar="character",   default='none',  help="Path to the folder containing the 10X folders"),
  make_option(c("-m", "--dataset_metadata_path"), type = "character",   metavar="character",   default='none',  help="Path to the Metadata matrix for each library (The first column should be named SampleID)"),
  make_option(c("-c", "--columns_metadata"),      type = "character",   metavar="character",   default='none',  help="Column names in the Metadata matrix (only factors allowed, not continuous variables)"),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output directory")
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
pkgs <- c("Seurat","dplyr","scales","RColorBrewer","biomaRt","ineq","vegan","rafalib")
inst_packages(pkgs)
#---------




### LOAD DATA AND SETUP Seurat OBJETC
#---------
cat("\nLoading/ data and metadata ...\n")
dataset_metadata <- as.data.frame(read.csv2(opt$dataset_metadata_path))
print(as.character(dataset_metadata[,1]))

datasets <- list.dirs(opt$input_path,recursive = F,full.names = F)
datasets <- datasets[datasets %in% as.character(dataset_metadata[,1])]

cat("\nThe following samples will be merged: ...\n")
print(datasets)

for(i in datasets ){
  cat("\nMerging sample: ",i," ...\n")
  a <- Read10X(paste0(opt$input_path,"/",i))
  colnames(a) <- paste0(colnames(a),"_",as.character(i))
  if(i != datasets[1]){
    temp <- CreateSeuratObject(a,project=i,min.cells = 0,is.expr = 0,min.genes = 0)
    DATA <- MergeSeurat(DATA, temp, do.normalize = FALSE)
  } else {
    DATA <- CreateSeuratObject(a,project=i,min.cells = 0,is.expr = 0,min.genes = 0 )
  }
}


#Add metadata
cat("\nThe following columns will be used ...\n")
use <- as.character(unlist(strsplit(opt$columns_metadata,","))) 
use <- use[use %in% colnames(dataset_metadata) ]
print(use)
for(i in use){
  DATA <- AddMetaData(object = DATA, metadata = setNames(dataset_metadata[match(DATA@meta.data$orig.ident, dataset_metadata[,1] ),i], rownames(DATA@meta.data)), col.name = i)}
#---------



### Saving the RAW Seurat object
#---------
cat("\nSaving the RAW Seurat object ...\n")
write.csv(DATA@meta.data,paste0(opt$output_path,"/QC_metadata_all_cells.csv"),row.names = T)
saveRDS(DATA, file = paste0(opt$output_path,"/Raw_Seurat_Object.rds") )
#---------



cat("\n!!! Script executed Sucessfully !!!\n")



### System and session information
#---------
cat("\n\n\n\n... SYSTEM INFORMATION ...\n")
Sys.info()

cat("\n\n\n\n... SESSION INFORMATION ...\n")
sessionInfo()
#---------