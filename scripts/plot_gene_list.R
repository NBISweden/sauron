#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
cat("\nRunning PLOTTING_GENE_LIST with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object"),
  make_option(c("-l", "--gene_list"),             type = "character",   metavar="character",   default='none',  help="Path to file contating the list of genes to be ploted. It should be a .csv file with one gene per column, one gene list per row."),
  make_option(c("-t", "--match_type"),            type = "character",   metavar="character",   default='exact',  help="The matching method to use. Use either 'exact' match or 'regex'. Defaut is exact."),
  make_option(c("-c", "--clustering_use"),        type = "character",   metavar="character",   default='none',  help="Column names in the Metadata matrix (only factors allowed, not continuous variables)"),
  make_option(c("-a", "--assay"),                 type = "character",   metavar="character",   default='RNA',  help="Assay to use for trajectory differential expression. Default is 'RNA'"),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output directory")
) 
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
#---------



##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")
initial.options <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(sub("--file=","",initial.options[grep("--file=",initial.options)]))
source( paste0(script_path,"/inst_packages.R") )
pkgs <- c("Seurat","dplyr","scales","RColorBrewer","rafalib")
inst_packages(pkgs)
#---------



#############################
### LOAD Seurat.v3 OBJECT ###
#############################
cat("\n### LOADING Seurat.v3 OBJECT ###\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------



##########################
### LOAD LIST OF GENES ###
##########################
cat("\nLoading/ list of genes ...\n")
gene_list <- read.csv2(opt$gene_list,header = T)
gene_list <- as.list(as.data.frame(gene_list))
gene_list <- lapply(gene_list, function(x) as.character(x[x!=""]) )
#---------


for(i in names(gene_list) ){
cat("\nProcessing list ",i," ...\n")
my_genes <- gene_list[[i]]
print(my_genes)

#########################################
### Grab genes present in the dataset ###
#########################################
if(casefold(opt$match_type) == 'exact'){
  sel <- casefold(my_genes) %in% casefold(rownames(DATA@assays[[opt$assay]]@data))
  if( sum(!sel) > 0 ){
    cat("\nThe following genes were not found in your dataset ... All other genes found will be plotted ...\n")
    print(my_genes[!sel])
  }
  my_genes <- unique(my_genes[sel])

} else {
  my_genes <- unique(rownames(DATA@assays[[opt$assay]]@data)[grepl(paste(my_genes,collapse = "|"),rownames(DATA@assays[[opt$assay]]@data))])
  if( length(my_genes) > 0 ){
    cat("\nThe following genes were found in your dataset using REGEX... \n")
    print(my_genes)
  }
}
#---------



############
### PLOT ###
############
if( length(my_genes) > 0 ){
  png(filename = paste0(opt$output_path,"/Heatmap_",i,".png"),width = 900,height = 900,res = 150)
  print(DoHeatmap(object = DATA, features = my_genes, assay=opt$assay, group.by = opt$clustering_use ))
  invisible(dev.off())
  
  png(filename = paste0(opt$output_path,"/ViolinPlot_",i,".png"),width = 200*3*10,height = 200*3*ceiling(length(as.character(unique(my_genes)))/10),res = 150)
  print(VlnPlot(object = DATA, features = as.character(my_genes),pt.size = .1,ncol=10,assay = opt$assay, group.by = opt$clustering_use))
  invisible(dev.off())
  
  for(j in names(DATA@reductions)){
    png(filename = paste0(opt$output_path,"/",j,"_plots_",i,".png"),width = 200*4*10,height = 200*3.5*ceiling(length(as.character(unique(my_genes)))/10),res = 150)
    print(FeaturePlot(object = DATA, features = as.character(my_genes), cols = c("grey", "blue"),reduction = j,pt.size = 1,ncol=10,dims = 1:2))
    invisible(dev.off())
  }
} else {
  cat("\nNone of the genes were found in your dataset... \n")
}
#---------
}



#############################
### SYSTEM & SESSION INFO ###
#############################
#---------
print_session_info()
#---------
