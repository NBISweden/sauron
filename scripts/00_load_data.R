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
cat("\nCREATING SEURAT OBJECT with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--input_path"),            type = "character",   metavar="character",   default='none',  help="Path to the folder containing the 10X folders"),
  make_option(c("-m", "--dataset_metadata_path"), type = "character",   metavar="character",   default='none',  help="Path to the Metadata matrix for each library (The first column should be named SampleID)"),
  make_option(c("-c", "--columns_metadata"),      type = "character",   metavar="character",   default='none',  help="Column names in the Metadata matrix (only factors allowed, not continuous variables)"),
  make_option(c("-c", "--integrate"),             type = "character",   metavar="character",   default='TRUE',  help="Logical specifying if dataset integration should be performed."),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output directory")
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
pkgs <- c("Seurat","dplyr","scales","RColorBrewer","biomaRt","ineq","vegan","rafalib")
inst_packages(pkgs)
#---------



#########################################
### LOAD DATA AND SETUP Seurat OBJETC ###
#########################################
cat("\nLoading/ data and metadata ...\n")
if(length(datasets) > 1){
  for(i in sort(datasets) ){
    a <- Read10X(paste0(opt$input_path,"/",i))
    colnames(a) <- paste0(colnames(a),"_",as.character(i))
    assign(i, CreateSeuratObject(a,project=i,min.cells = 0,min.features = 0))
  }
  DATA <- merge(sort(datasets)[1],y=sort(datasets)[-1])
} else {
  a <- Read10X(paste0(opt$input_path,"/",i))
  colnames(a) <- paste0(colnames(a),"_",as.character(i))
  DATA <- CreateSeuratObject(a,project=i,min.cells = 0,min.features = 0)
}
rm(c(a,datasets)); gc()
#---------



####################
### ADD METADATA ###
####################
cat("\nThe following columns will be used ...\n")
use <- as.character(unlist(strsplit(opt$columns_metadata,","))) 
use <- use[use %in% colnames(dataset_metadata) ]
print(use)
for(i in use){
  DATA <- AddMetaData(object = DATA, metadata = setNames(dataset_metadata[match(DATA@meta.data$orig.ident, dataset_metadata[,1] ),i], rownames(DATA@meta.data)), col.name = i)}
#---------



####################################
### INTEGRATE DATASETS USING CCA ###
####################################
cat("\nIntegrating datasets with CCA ...\n")
if( as.logical(opt$integrate) & (length(datasets) > 1) ){
  DATA.list <- SplitObject(DATA, split.by = "orig.ident")
  for (i in 1:length(DATA.list)) {
    DATA.list[[i]] <- NormalizeData(DATA.list[[i]], verbose = FALSE)
    DATA.list[[i]] <- FindVariableFeatures(DATA.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }
  DATA.anchors <- FindIntegrationAnchors(object.list = DATA.list, dims = 1:30)
  DATA.integrated <- IntegrateData(anchorset = DATA.anchors, dims = 1:30)
  DefaultAssay(DATA.integrated) <- "integrated"
}
#---------



####################################
### Saving the RAW Seurat object ###
####################################
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