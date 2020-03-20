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
  make_option(c("-t", "--estimate_molecules_from_read_count"),     type = "character",   metavar="character",   default='none',  help="Whether to estimate mRNA molecules from the read counts. This makes the data more comparable to UMI counts from Drop-seq."),
  make_option(c("-n", "--mapping_threshold"),     type = "character",   metavar="character",   default='2.5',  help="The lower limit threshold for counting a molecule in range 1-100. 1 equals at least 1 read per gene. Default is 2.5 reads per gene."),
  make_option(c("-g", "--sum_to_gene_level"),     type = "character",   metavar="character",   default='none',  help="Whether to summarize Ensenbl gene IDS to gene level"),
  make_option(c("-s", "--species_use"),           type = "character",   metavar="character",   default='hsapiens',  help="Ensembl species to use"),
  make_option(c("-a", "--assay"),                 type = "character",   metavar="character",   default='RNA',   help="Assay to be used in the analysis."),
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
#pkgs <- c("Seurat","dplyr","rafalib","Matrix","parallel")
#inst_packages(pkgs)

suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(rafalib)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(parallel)))
suppressMessages(suppressWarnings(library(biomaRt)))

#---------



############################################
### LOAD DATA AND SETUP Seurat.v3 OBJECT ###
############################################
cat("\nLoading/ data and metadata ...\n")

L <- readLines(opt$dataset_metadata_path, n=1)
if (grepl(";", L)) {
  dataset_metadata <- as.data.frame(read.csv2(opt$dataset_metadata_path))
} else {
  dataset_metadata <- as.data.frame(read.csv(opt$dataset_metadata_path))
}
print(as.character(dataset_metadata[,1]))

datasets <- list.dirs(opt$input_path,recursive = F,full.names = F)
datasets <- sort(datasets[datasets %in% as.character(dataset_metadata[,1])])

cat("\nThe following samples will be merged: ...\n")
print(datasets)

if(length(datasets) > 1){
  #for(i in sort(datasets) ){
  cat("\nloading datasets\n")
  cl <- makeCluster(detectCores()-1,type = "FORK")
  clusterExport(cl, varlist = c("datasets","opt") )
  data <- parLapplyLB(cl, datasets, function(i){
  # data <- lapply(datasets, function(i){
    cat("Processing dataset ",i)
    require(Seurat)
    require(Matrix)
    require(utils)
    if( sum(grepl(".mtx", list.files(paste0(opt$input_path,"/",i)))) >= 1 ){
      #read 10X files
      a <- NULL
      try(a <- Seurat::Read10X(paste0(opt$input_path,"/",i)),silent = T)
      print(dim(a))
      if(is.null(a)){
        a <- as(Matrix::readMM(paste0(opt$input_path,"/",i,"/matrix.mtx")),Class = "dgCMatrix")
        colnames(a) <- as.character(read.table(paste0(opt$input_path,"/",i,"/barcodes.tsv"))[,1])
        rownames(a) <- as.character(read.table(paste0(opt$input_path,"/",i,"/features.tsv"))[,1])
      }
      
    } else if  ( sum(grepl(".h5", list.files(paste0(opt$input_path,"/",i)))) == 1 ){

      a <- Seurat::Read10X_h5( paste0(opt$input_path,"/",i,"/", grep(".h5", list.files(paste0(opt$input_path,"/",i)),value = T)  ) )
      if(is.null(dim(a))){a <- a[[1]]}

    } else if  ( sum(grepl(".csv", list.files(paste0(opt$input_path,"/",i)))) == 1 ){
      #read .csv files
      a <- read.csv2(paste0(opt$input_path,"/",i,"/",grep(".csv", list.files(paste0(opt$input_path,"/",i)),value = T) ),row.names = 1 )
      if(ncol(a) == 0){a <- read.csv(paste0(opt$input_path,"/",i,"/",grep(".csv", list.files(paste0(opt$input_path,"/",i)),value = T) ),row.names = 1 )}
      a <- Matrix::Matrix(as.matrix(rowsum(a,sub("[_.,].*","",rownames(a)))),sparse=T)
    
    } else if  ( sum(grepl(".tsv|.txt", list.files(paste0(opt$input_path,"/",i)))) == 1 ){
      #read .tsv / .txt files
      a <- read.delim(paste0(opt$input_path,"/",grepl(".txt|.tsv", list.files(paste0(opt$input_path,"/",i)),value = T) ),row.names = T ,header = T)
      a <- Matrix::Matrix(as.matrix(rowsum(a,sub("[_.,].*","",rownames(a)))),sparse=T)
    }
    
    colnames(a) <- paste0(sub("-.*","",colnames(a)),"_",as.character(i))
    #assign(i, CreateSeuratObject(a,project=i,min.cells = 1,min.features = 1),envir = .GlobalEnv)
    cat("The size of dataset", i, " is: ", dim(a),"\n" )
    return(a)
  })
  names(data) <- datasets
  cat("\nDimension of loaded datasets\n")
  print(as.data.frame(lapply(data,dim),row.names = c("genes","cells")))

  #}
  
  cat("Merging datasets\n" )
  #DATA <- merge(get(sort(datasets)[1]), y=mget(sort(datasets)[-1]),all=T)
  all_genes <- unique(unlist(parLapplyLB(cl,data,function(x) return(rownames(x)))))

  clusterExport(cl, varlist = c("all_genes") )
  data <- parLapplyLB(cl, data, all_genes=all_genes,function(x,all_genes) {
    m <- Matrix::Matrix(0,nrow = length(all_genes), ncol = ncol(x),sparse = T,dimnames = list(all_genes,colnames(x)))
    m[rownames(x),] <- x
    return(m)
  })
  
  cat("New dimensions\n" )
  print(as.data.frame(lapply(data,dim),row.names = c("genes","cells")))
  DATA <- do.call(cbind,data)
  
  mart = useMart("ensembl", dataset = paste0(opt$species_use,"_gene_ensembl"),host="jul2019.archive.ensembl.org")
  annot <- getBM(c("ensembl_gene_id","external_gene_name","transcript_length","gene_biotype","chromosome_name"),mart = mart)
  
  if(casefold(opt$estimate_molecules_from_read_count) %in% c("yes","true") ){
    cat("Calculating estimate_molecules\n" )
    item <- annot[match(rownames(DATA) , annot[,1]),3]
    DATA <- DATA[!is.na(item),] / item[!is.na(item)]
    DATA <- round( DATA * 100 / as.numeric(opt$mapping_threshold) )
  }
  
  
  if(casefold(opt$sum_to_gene_level) %in% c("yes","true") ){
    cat("Summing up to gene level\n" )
    item <- annot[match(rownames(DATA) , annot[,1]),2]
    item[is.na(item)] <- "NA"
    item[item==""] <- "NA"
    item[grep("-AS1$",item)] <- "NA"
    DATA <- Matrix::Matrix( rowsum( as.matrix(DATA), item ),sparse = T)
    DATA <- DATA[ rownames(DATA) != "NA" , ]
  }
  
  cat("Creating Seurat Object\n" )
  DATA <- CreateSeuratObject(DATA,min.cells = 1,min.features = 1,assay = opt$assay)
  DATA$orig.ident <- setNames(sub("(.*?)_","",colnames(DATA)) , colnames(DATA) )
  rm(data); invisible(gc())
  
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
use <- colnames(dataset_metadata)
print(use)
for(i in use){
  DATA <- AddMetaData(object = DATA, 
                      metadata = setNames(as.character(dataset_metadata[match(as.character(DATA$orig.ident), as.character(dataset_metadata[,1]) ),i]), rownames(DATA@meta.data)), col.name = i)}


#Load individual metadata files from each sample
###
cat("Adding individual individual metadata files from each sample\n" )
# indiv_metadata <- parLapplyLB(cl, datasets, function(i){
indiv_metadata <- lapply(datasets, function(i){
  #data <- lapply(datasets, function(i){
  cat("\nProcessing metadata from dataset ",i)
  if( sum(grepl("metadata", list.files(paste0(opt$input_path,"/",i)))) == 1 ){
    a <- read.csv2(paste0(opt$input_path,"/",i,"/",grep("metadata", list.files(paste0(opt$input_path,"/",i)),value = T) ),row.names = 1 )
    if(ncol(a) == 0){a <- read.csv(paste0(opt$input_path,"/",i,"/",grep("metadata", list.files(paste0(opt$input_path,"/",i)),value = T) ),row.names = 1 )}
    
    rownames(a) <- paste0(rownames(a),"_",as.character(i))
    return(a)
  }
})


cat("Merging metadata\n" )
cat("The following unique individual metadata variables were found across all datasets:\n" )
print(unique(unlist(lapply(indiv_metadata,colnames))))

#DATA <- merge(get(sort(datasets)[1]), y=mget(sort(datasets)[-1]),all=T)
all_metadata <- unique( c( colnames(DATA@meta.data), unlist(lapply(indiv_metadata,colnames)) ) )
dat <- data.frame(matrix( NA , nrow = ncol(DATA), ncol = length(all_metadata), dimnames = list(colnames(DATA), all_metadata)),stringsAsFactors = F)
# dat <- droplevels(dat)
# dat[ rownames(DATA@meta.data), colnames(DATA@meta.data) ] <- droplevels(DATA@meta.data)

cat("Adding individual metadata\n" )

for(j in 1:length(indiv_metadata) ){
  temp <- indiv_metadata[[j]]
  common <- rownames(temp) [ rownames(temp) %in% rownames(dat)  ]
  for(i in colnames(temp)){
    dat[  common ,i] <- as.character( temp[common,i] ) }
}

cat("Adding overall metadata\n" )
for(i in colnames(DATA@meta.data)){
  dat[rownames(DATA@meta.data),i] <- as.character(DATA@meta.data[,i]) }

DATA@meta.data <- dat
rm(dat)
#---------



###################################
### SAVING RAW Seurat.v3 OBJECT ###
###################################
cat("\nSaving the RAW Seurat object ...\n")
write.csv2(DATA@meta.data,paste0(opt$output_path,"/QC_metadata_all_cells.csv"),row.names = T)
saveRDS(DATA, file = paste0(opt$output_path,"/raw_seurat_object.rds") )
stopCluster(cl)
#---------



#############################
### SYSTEM & SESSION INFO ###
#############################
#---------
print_session_info()
#---------




