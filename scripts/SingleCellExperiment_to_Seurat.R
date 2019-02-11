#!/usr/bin/env Rscript

### LOAD LIBRARIES
#---------
library(optparse)
#---------


### DEFINE PATH TO LOCAL FILES
#---------
cat("\nRunning SingleCellExperiment to Seurat conversion with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--input_file"),       type = "character",   metavar="character",   default='none',  help="Path to RDS file contating the SingleCellExperiment object"),
  make_option(c("-s", "--species_use"),           type = "character",   metavar="character",   default='none',  help="Species from the sample for cell scoring"),
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
pkgs <- c("Seurat","dplyr","scales","RColorBrewer","rafalib","biomaRt","SingleCellExperiment","stats4")
inst_packages(pkgs)
#---------




### LOAD SingleCellExperiment OBJECT 
#---------
cat("\nLoading  ...\n")
DATA <- readRDS(opt$input_file)
cat("\nData loaded! Extracting information from the object ...\n")

metadata <- colData(DATA)
data <- as.matrix(counts(DATA))

cat("\nGathering information on ensembl gene IDs \n")
if(opt$species_use == "mouse"){
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  annot <- getBM(c("ensembl_gene_id","mgi_symbol"), mart=mouse)
  a <- annot[match(rownames(data),annot[,"ensembl_gene_id"]),"mgi_symbol"]
} else if (opt$species_use == "human"){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  annot <- getBM(c("ensembl_gene_id","hgnc_symbol"), mart=human)
  a <- annot[match(rownames(data),annot[,"ensembl_gene_id"]),"hgnc_symbol"]
}
#---------



#Annotate and sum counts per gene
#---------
cat("\nAnnotating ensembl gene IDs to gene names ...\n")
a[is.na(a)] <- ""
data <- rowsum(data, group = a)
data <- data[-1,]

cat("\nCreating the RAW Seurat object ...\n")
DATA <- CreateSeuratObject(raw.data = data,min.cells = 0,min.genes = 0,is.expr = 0,meta.data = as.data.frame(metadata) )
#---------





### Saving the RAW Seurat object
#---------
cat("\nSaving the RAW Seurat object ...\n")
write.csv(DATA@meta.data,paste0(opt$output_path,"/QC_metadata_all_cells.csv"),row.names = T)
saveRDS(DATA, file = paste0(opt$output_path,"/SCExp_Seurat_Object.rds") )
#---------



cat("\n!!! Script executed Sucessfully !!!\n")



### System and session information
#---------
cat("\n\n\n\n... SYSTEM INFORMATION ...\n")
Sys.info()

cat("\n\n\n\n... SESSION INFORMATION ...\n")
sessionInfo()
#---------