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
  make_option(c("-b", "--integrate"),             type = "character",   metavar="character",   default='TRUE',  help="Logical specifying if dataset integration should be performed."),
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
pkgs <- c("Seurat","dplyr","rafalib","Matrix","parallel")
inst_packages(pkgs)
#---------



############################################
### LOAD DATA AND SETUP Seurat.v3 OBJECT ###
############################################
cat("\nLoading/ data and metadata ...\n")
dataset_metadata <- as.data.frame(read.csv2(opt$dataset_metadata_path))
print(as.character(dataset_metadata[,1]))

datasets <- list.dirs(opt$input_path,recursive = F,full.names = F)
datasets <- sort(datasets[datasets %in% as.character(dataset_metadata[,1])])

cat("\nThe following samples will be merged: ...\n")
print(datasets)

if(length(datasets) > 1){
  #for(i in sort(datasets) ){
  cat("\nloading datasets\n")
  cl <- makeCluster(detectCores()-1,type = "FORK")
  clusterEvalQ(cl, library(rms))
  clusterExport(cl, varlist = c("datasets","opt") )
  data <- parLapplyLB(cl, datasets, function(i){
    require(Seurat)
    require(Matrix)
    require(utils)
    if( sum(grepl(".mtx", list.files(paste0(opt$input_path,"/",i)))) == 1 ){
      #read 10X files
      a <- Seurat::Read10X(paste0(opt$input_path,"/",i))
    } else if  ( sum(grepl(".csv", list.files(paste0(opt$input_path,"/",i)))) == 1 ){
      #read .csv files
      a <- read.csv(paste0(opt$input_path,"/",i,"/",grep(".csv", list.files(paste0(opt$input_path,"/",i)),value = T) ),row.names = 1 )
      if(ncol(a) == 0){a <- read.csv2(paste0(opt$input_path,"/",i,"/",grep(".csv", list.files(paste0(opt$input_path,"/",i)),value = T) ),row.names = 1 )}
      a <- Matrix::Matrix(as.matrix(rowsum(a,sub("[_.,].*","",rownames(a)))),sparse=T)
      } else if  ( sum(grepl(".tsv|.txt", list.files(paste0(opt$input_path,"/",i)))) == 1 ){
      #read .csv files
      a <- read.delim(paste0(opt$input_path,"/",grepl(".txt|.tsv", list.files(paste0(opt$input_path,"/",i)),value = T) ),row.names = T ,header = T)
      a <- Matrix::Matrix(as.matrix(rowsum(a,sub("[_.,].*","",rownames(a)))),sparse=T)
      }
    
    colnames(a) <- paste0(colnames(a),"_",as.character(i))
    #assign(i, CreateSeuratObject(a,project=i,min.cells = 1,min.features = 1),envir = .GlobalEnv)
    cat("The size of dataset", i, " is: ", dim(a),"\n" )
    return(a)
  })
  
  cat("\nDimension of loaded datasets\n")
  print(lapply(data,dim))
  
  #}
  cat("Merging datasets\n" )
  #DATA <- merge(get(sort(datasets)[1]), y=mget(sort(datasets)[-1]),all=T)
  all_genes <- unique(unlist(parLapplyLB(cl,data,function(x) return(rownames(x)))))
  print(length(all_genes))
  
  clusterExport(cl, varlist = c("all_genes") )
  data <- parLapplyLB(cl, data, all_genes=all_genes,function(x,all_genes) {
    m <- Matrix::Matrix(0,nrow = length(all_genes), ncol = ncol(x),sparse = T,dimnames = list(all_genes,colnames(x)))
    m[rownames(x),] <- x
    return(m)
  })
  
  DATA <- do.call(cbind,data)
  rm(data); invisible(gc())
  DATA <- CreateSeuratObject(DATA,min.cells = 1,min.features = 1)
} else {
  a <- Read10X(paste0(opt$input_path,"/",i))
  colnames(a) <- paste0(colnames(a),"_",as.character(i))
  DATA <- CreateSeuratObject(a,min.cells = 1,min.features = 1)
}
cat("\nThe total dimensions of your merged raw dataset is: ",dim(DATA),"\n")
invisible(gc())
#---------



####################
### ADD METADATA ###
####################
cat("\nThe following columns will be used ...\n")
use <- as.character(unlist(strsplit(opt$columns_metadata,","))) 
use <- use[use %in% colnames(dataset_metadata) ]
print(use)
for(i in use){
  DATA <- AddMetaData(object = DATA, metadata = setNames(dataset_metadata[match(DATA$orig.ident, dataset_metadata[,1] ),i], rownames(DATA@meta.data)), col.name = i)}
#---------



###################################
### SAVING RAW Seurat.v3 OBJECT ###
###################################
cat("\nSaving the RAW Seurat object ...\n")
write.csv(DATA@meta.data,paste0(opt$output_path,"/QC_metadata_all_cells.csv"),row.names = T)
saveRDS(DATA, file = paste0(opt$output_path,"/Raw_Seurat_Object.rds") )
stopCluster(cl)
#---------



#############################
### SYSTEM & SESSION INFO ###
#############################
#---------
print_session_info()
#---------




