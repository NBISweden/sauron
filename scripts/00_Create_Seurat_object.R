#!/usr/bin/env Rscript

### LOAD LIBRARIES
#---------
library(optparse)
#---------


### DEFINE PATH TO LOCAL FILES
#---------
cat("\nRunning scQC with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--TenX_files_path"),       type = "character",   metavar="character",   default='none',  help="Path to the folder containing the 10X folders"),
  make_option(c("-m", "--dataset_metadata_path"), type = "character",   metavar="character",   default='none',  help="Path to the Metadata matrix for each library (The first column should be named SampleID)"),
  make_option(c("-c", "--columns_metadata"),      type = "character",   metavar="character",   default='none',  help="Column names in the Metadata matrix (only factors allowed, not continuous variables)"),
  make_option(c("-f", "--aux_functions_path"),    type = "character",   metavar="character",   default='none',  help="File with supplementary functions"),
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
source(opt$aux_functions_path)
pkgs <- c("Seurat","dplyr","scales","RColorBrewer","biomaRt","ineq","vegan","rafalib")
inst_packages(pkgs)
#---------




### LOAD DATA AND SETUP Seurat OBJETC
#---------
cat("\nLoading/ data and metadata ...\n")
dataset_metadata <- as.data.frame(read.csv2(opt$dataset_metadata_path))
head(dataset_metadata)

datasets <- list.dirs(opt$TenX_files_path)[-1]
print(datasets)

for(i in 1:length(datasets) ){
  a <- Read10X(datasets[i])
  colnames(a) <- paste0(colnames(a),"_",as.character(dataset_metadata$SampleID[i]))
  if(i>1){
    temp <- CreateSeuratObject(a,project=as.character(dataset_metadata$SampleID[i]),min.cells = 10,is.expr = 1,min.genes = 200)
    DATA <- MergeSeurat(DATA, temp, do.normalize = FALSE)
  } else {
    DATA <- CreateSeuratObject(a,project=as.character(dataset_metadata$SampleID[i]),min.cells = 10,is.expr = 1,min.genes = 200 )
  }
}

#Add metadata
cat("\nThe following columns will be used ...\n")
use <- as.character(unlist(strsplit(opt$columns_metadata,",")))
print(use)
for(i in use){
DATA <- AddMetaData(object = DATA, metadata = setNames(dataset_metadata[match(DATA@meta.data$orig.ident, dataset_metadata$SampleID),i], rownames(DATA@meta.data)), col.name = i)}
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