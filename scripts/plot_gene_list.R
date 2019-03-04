#!/usr/bin/env Rscript

### LOAD LIBRARIES
#---------
library(optparse)
#---------


### DEFINE PATH TO LOCAL FILES
#---------
cat("\nRunning scQC with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object"),
  make_option(c("-l", "--gene_list"),             type = "character",   metavar="character",   default='none',  help="Path to file contating the list of genes to be ploted. It should be a .csv file with one gene per column, one gene list per row."),
  make_option(c("-t", "--match_type"),            type = "character",   metavar="character",   default='exact',  help="The matching method to use. Use either 'exact' match or 'regex'. Defaut is exact."),
  make_option(c("-c", "--clustering_use"),      type = "character",   metavar="character",   default='none',  help="Column names in the Metadata matrix (only factors allowed, not continuous variables)"),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output directory")
) 
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
#---------




### LOAD LIBRARIES
#---------
cat("\nLoading/installing libraries ...\n")
initial.options <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(sub("--file=","",initial.options[grep("--file=",initial.options)]))
source( paste0(script_path,"/inst_packages.R") )
pkgs <- c("Seurat","dplyr","scales","RColorBrewer","rafalib")
inst_packages(pkgs)
#---------


### LOAD Seurat OBJECT 
#---------
cat("\nLoading/ data and metadata ...\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------




### LOAD list of genes
#---------
cat("\nLoading/ list of genes ...\n")
gene_list <- read.csv2(opt$gene_list,header = F)
rownames(gene_list) <- as.character(gene_list[,1])
gene_list <- as.list(as.data.frame(t(gene_list[,-1])))
gene_list <- lapply(gene_list, function(x) as.character(x[x!=""]) )
#---------


for(i in names(gene_list) ){
cat("\nProcessing list ",i," ...\n")
  
my_genes <- gene_list[[i]]
print(my_genes)

### Grab genes present in the dataset
#---------
if(casefold(opt$match_type) == 'exact'){
  sel <- my_genes %in% rownames(DATA@scale.data)
  if( sum(!sel) > 0 ){
    cat("\nThe following genes were not found in your dataset ... All other genes found will be plotted ...\n")
    print(my_genes[!sel])
  }
  my_genes <- unique(my_genes[sel])

} else {
  my_genes <- unique(rownames(DATA@scale.data)[grepl(paste(my_genes,collapse = "|"),rownames(DATA@scale.data))])
  if( length(my_genes) > 0 ){
    cat("\nThe following genes were found in your dataset using REGEX... \n")
    print(my_genes)
  }
}
#---------



### Plot 
#---------
if( length(my_genes) > 0 ){
  if( opt$clustering_use != "none" ){  
    DATA@ident <- factor(NULL)
    DATA <- SetIdent(DATA,ident.use = DATA@meta.data[,opt$clustering_use])
  }
  
  png(filename = paste0(opt$output_path,"/Heatmap_",i,".png"),width = 900,height = 900,res = 150)
  print(DoHeatmap(object = DATA, genes.use = my_genes, slim.col.label = TRUE, remove.key = TRUE, use.scaled = F))
  invisible(dev.off())
  
  png(filename = paste0(opt$output_path,"/ViolinPlot_",i,".png"),width = 200*3*10,height = 200*3*ceiling(length(as.character(unique(my_genes)))/10),res = 150)
  print(VlnPlot(object = DATA, features.plot = as.character(my_genes),point.size.use = .1,nCol=10))
  invisible(dev.off())
  
  png(filename = paste0(opt$output_path,"/tSNE_plots_",i,".png"),width = 200*3*10,height = 200*3*ceiling(length(as.character(unique(my_genes)))/10),res = 150)
  print(FeaturePlot(object = DATA, features.plot = as.character(my_genes), cols.use = c("grey", "blue"), reduction.use = "tsne",pt.size = 1,nCol=10))
  invisible(dev.off())
  
} else {
  cat("\nNone of the genes were found in your dataset... \n")
}
#---------
}


cat("\n!!! Script executed Sucessfully !!!\n")


### System and session information
#---------
cat("\n\n\n\n... SYSTEM INFORMATION ...\n")
INFORMATION <- Sys.info()
print(as.data.frame(INFORMATION))

cat("\n\n\n\n... SESSION INFORMATION ...\n")
sessionInfo()
#---------

