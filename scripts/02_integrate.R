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
cat("\nRunning DATA INTEGRATION with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object"),
  make_option(c("-c", "--columns_metadata"),      type = "character",   metavar="character",   default='none',  help="Column names in the Metadata matrix (only factors allowed, not continuous variables)"),
  make_option(c("-r", "--regress"),               type = "character",   metavar="character",   default='none',  help="Variables to be regressed out using linear modeling."),
  make_option(c("-b", "--integration_method"),    type = "character",   metavar="character",   default='cca,orig.ident',  help="Integration method to be used. 'cca', mmn', 'scale' and 'combat' are available at the moment. The batches (column names in the metadata matrix) to be removed should be provided as arguments comma separated. E.g.: 'Combat,sampling_day'. For MNN, an additional integer parameter is supplied as the k-nearest neighbour."),
  make_option(c("-v", "--var_genes"),             type = "character",   metavar="character",   default='scran,.2',  help="Whether use 'Seurat' or the 'Scran' method for variable genes identification. An additional value can be placed after a comma to define the level of dispersion wanted for variable gene selection. 'Seurat,2' will use the threshold 2 for gene dispersions. Defult is 'scran,0.001'. For Scran, the user should inpup the level of biological variance 'Scran,0.001'. An additional blocking parameter (a column from the metadata) can ba supplied to 'Scran' method block variation comming from uninteresting factors, which can be parsed as 'Scran,0.2,Batch'."),
  make_option(c("-s", "--cluster_use"),           type = "character",   metavar="character",   default='all',  help="The cluster to be used for analysis."),
  make_option(c("-a", "--assay"),                 type = "character",   metavar="character",   default='RNA',  help="The default assay to use to integrate."),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output directory")
) 
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)

col_scale <- c("grey85","navy")
VAR_choice <- as.character(unlist(strsplit(opt$var_genes,",")))
#---------



###################################
### SETUP MULTICORE ENVIRONMENT ###
###################################
#options(future.globals.maxSize= 2048 * 1024^2)
#plan(strategy = "multicore", workers = nbrOfWorkers())
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
#inst_packages(pkgs)

suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(scales)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
suppressMessages(suppressWarnings(library(biomaRt)))
suppressMessages(suppressWarnings(library(igraph)))
suppressMessages(suppressWarnings(library(sva)))
suppressMessages(suppressWarnings(library(rafalib)))
suppressMessages(suppressWarnings(library(parallel)))
suppressMessages(suppressWarnings(library(scran)))
suppressMessages(suppressWarnings(library(scater)))
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
  cells_use <- rownames(DATA@meta.data)}
# sel <- rowSums(as.matrix(DATA@assays[[opt$assay]]@counts) >= 1) >= 1
# DATA <- CreateSeuratObject(as.matrix(DATA@assays[[opt$assay]]@counts[sel,cells_use]), meta.data = DATA@meta.data[cells_use,])
#---------



######################
### NORMALIZE DATA ###
######################
#cat("\nNormalizing and identifying highly variable genes ...\n")
#DATA <- NormalizeData(DATA,scale.factor = 1000)
#---------



#######################################
### INTEGRATE DATASETS USING COMBAT ###
#######################################
integration_method <- unlist(strsplit(opt$integration_method,","))

if ((length(integration_method) >= 2) & (casefold(integration_method[1]) == "combat") ){
  cat("\n### INTEGRATING DATASETS USING COMBAT ###\n")
  print(integration_method)
  
  #Defining batch variables
  batch <- factor(DATA@meta.data[,integration_method[2]])
  mod0 <- model.matrix(~1, data=as.data.frame(DATA@meta.data))
  
  #Transforming counts to log
  logdata <- log2(as.matrix(DATA@assays[[opt$assay]]@data)[,rownames(DATA@meta.data)]+1)
  sum(rowSums(logdata) == 0)
  logdata <- logdata[rowSums(logdata) != 0,]
  
  #Applying Combat
  combat_data <- ComBat(dat=logdata, batch=batch, mod=mod0)
  
  #Transforming counts back to original scale (and removing negative values to 0)
  combat_data <- round(2^(combat_data)-1,0)
  sum(combat_data < 0)
  combat_data[combat_data < 0] <- 0
  DATA@assays[["combat"]] <- CreateAssayObject(data = combat_data,min.cells = 0,min.features = 0)
  DATA <- NormalizeData(DATA,scale.factor = 1000)
  rm(combat_data,logdata,mod0);  invisible(gc())
}
# if( prod(dim(DATA@assays[[opt$assay]]@data) == c(0,0))!=0 ){ DATA <- var_gene_method(DATA,VAR_choice) }
#---------



####################################
### INTEGRATE DATASETS USING MNN ###
####################################
if ((length(integration_method) >= 2) & (casefold(integration_method[1]) == "mnn") ){
  cat("\n### INTEGRATING DATASETS USING MNN ###\n")
  print(integration_method)
  
  #Defining batch variables
  cat("\nCreatting dataset list\n")
  batch <- as.character(factor(DATA@meta.data[,integration_method[2]]))
  
  DATA.list <- SplitObject(DATA, split.by = integration_method[2])
  
  #DATA.list <- lapply(DATA.list, function(x){ 
  #  })

  if( (length(DATA.list) > 1) ){
    
    # define HVGs per dataset
    cat("\nComputing HVGs\n")
    for (i in 1:length(DATA.list)) {
      cat("\nProcessing dataset: ",names(DATA.list)[i]," \n")
      #DATA.list[[i]] <- NormalizeData(DATA.list[[i]], verbose = FALSE,scale.factor = 1000)
      if(file.exists(paste0(opt$output_path,"/variable_genes/var_genes_",names(DATA.list)[i],"/HVG_info_",VAR_choice[1],".csv"))){
        cat("\nVariable genes for this dataset found, and will be used. Skiping HVG calculation.\n")
        temp <- read.csv2(paste0(opt$output_path,"/variable_genes/var_genes_",names(DATA.list)[i],"/HVG_info_",VAR_choice[1],".csv"),row.names=1)
        DATA.list[[i]]@assays[[opt$assay]]@meta.features <- temp
        DATA.list[[i]]@assays[[opt$assay]]@var.features <- rownames(temp)[temp$use]
      } else {DATA.list[[i]] <- compute_hvgs(DATA.list[[i]],VAR_choice,paste0(opt$output_path,"/variable_genes/var_genes_",names(DATA.list)[i]),assay = opt$assay)}
    }

    # select the most informative genes that are shared across all datasets:
    #universe <- Reduce(intersect, lapply(DATA.list,function(x){x@assays[[opt$assay]]@var.features}))
    cat("\nComputing HVGs\n")
    universe <- unique(unlist(lapply(DATA.list,function(x){x@assays[[opt$assay]]@var.features})))
    head(universe,50)
    cat("\n",length(universe)," genes found as variable within datasets\n")
    
    #Separating batch matricies
    DATA.list <- lapply(DATA.list, function(x){ x@assays[[opt$assay]]@data[universe,]} )
    #rm(DATA.list)
    #datasets <- names(myinput)
    myinput <- list()
    
    if(  is.na(integration_method[3]) ) { myinput[["k"]] <- 20
    } else { myinput[["k"]] <- as.numeric(integration_method[3]) }
    myinput[["approximate"]] <-  TRUE
    myinput[["d"]] <-  51
    #myinput[["BPPARAM"]] <-  MulticoreParam(workers = detectCores()-1)

    #Applying MNN correction on raw counts
    cat("\nApplying MNN correction on normalized counts\n")
    out <- do.call(fastMNN,args = c(DATA.list,myinput) )
    cat("\nMNN computation done\n")
    out <- t(out$corrected)
    colnames(out) <- unlist(lapply(DATA.list,function(x){colnames(x)}))
    out <- out[,colnames(DATA)]
    rownames(out) <- paste0("dim",1:myinput$d)
    #DATA@assays[["mnn"]] <- CreateAssayObject(data = out,min.cells = 0,min.features = 0)
    DATA@reductions[["mnn"]] <- CreateDimReducObject(embeddings = t(out),key = "MNN_",assay = opt$assay)
    #DATA@assays$mnn@var.features <- rownames(DATA@assays$mnn@data)
    rm(out, myinput);  invisible(gc())
  }
}
#---------





####################################
### INTEGRATE DATASETS USING CCA ###
####################################
if ((length(integration_method) >= 1) & (casefold(integration_method[1]) == "cca") ){
  cat("\n### INTEGRATING DATASETS USING CCA ###\n")
  print(integration_method)
  
  DATA.list <- SplitObject(DATA, split.by = integration_method[2])
  if( (length(DATA.list) > 1) ){
    
    # DATA.list <- lapply(DATA.list,function(x){
    #   #x <- NormalizeData(x, verbose = FALSE,scale.factor = 1000)
    #   x <- compute_hvgs(x,VAR_choice,paste0(opt$output_path,"/var_genes_",names(DATA.list)[x]),assay = opt$assay)
    #   return(x)
    # })
    # 
    for (i in 1:length(DATA.list)) {
      #DATA.list[[i]] <- NormalizeData(DATA.list[[i]], verbose = FALSE,scale.factor = 1000)
      DATA.list[[i]] <- compute_hvgs(DATA.list[[i]],VAR_choice,paste0(opt$output_path,"/var_genes_",names(DATA.list)[i]),assay = opt$assay)
      gc()
    }
    universe <- unique(unlist(lapply(DATA.list,function(x){x@assays[[opt$assay]]@var.features})))

    DATA.anchors <- FindIntegrationAnchors(object.list = DATA.list, dims = 1:30, anchor.features = universe)
    DATA <- IntegrateData(anchorset = DATA.anchors, dims = 1:30, new.assay.name = "cca")
    DATA@assays$cca@var.features <- rownames(DATA@assays$cca@data)
    rm(DATA.list,DATA.anchors); gc()
  }
}
#---------



###########################
### FIND VARIABLE GENES ###
###########################
if(DefaultAssay(DATA) == opt$assay){
  output_path <- paste0(opt$output_path,"/variable_genes/All_datasets_together")
  DATA <- compute_hvgs(DATA,VAR_choice,output_path,assay = opt$assay)}
#---------



###################################
### SAVING Seurat.v3 OBJECT ###
###################################
cat("\n### Saving Seurat object ###\n")
saveRDS(DATA, file = paste0(opt$output_path,"/seurat_object.rds") )
#---------



#############################
### SYSTEM & SESSION INFO ###
#############################
print_session_info()
#---------

