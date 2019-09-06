#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
cat("\nRunning QUALITY CONTROL with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object"),
  make_option(c("-m", "--columns_metadata"),      type = "character",   metavar="character",   default='none',  help="Column names in the Metadata matrix (only factors allowed, not continuous variables)"),
  make_option(c("-s", "--species_use"),           type = "character",   metavar="character",   default='none',  help="Species from the sample for cell scoring"),
  make_option(c("-n", "--remove_non_coding"),     type = "character",   metavar="character",   default='True',     help="Removes all non-coding and pseudogenes from the data. Default is 'True'."),
  make_option(c("-p", "--plot_gene_family"),      type = "character",   metavar="character",   default='Rps,Rpl,mt-,Hb',  help="Gene families to plot QC. They should start with the pattern."),
  make_option(c("-r", "--remove_gene_family"),    type = "character",   metavar="character",   default='mt-',  help="Gene families to remove from the data after QC. They should start with the pattern."),
  make_option(c("-g", "--min_gene_count"),        type = "character",   metavar="character",   default='5',  help="Minimun number of cells needed to consider a gene as expressed. Defaults to 5."),
  make_option(c("-c", "--min_gene_per_cell"),        type = "character",   metavar="character",   default='200', help="Minimun number of genes in a cell needed to consider a cell as good quality. Defoust to 200."),
  make_option(c("-a", "--assay"),                 type = "character",   metavar="character",   default='rna',   help="Assay to be used in the analysis."),
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
pkgs <- c("Seurat","dplyr","scales","RColorBrewer","biomaRt","ineq","vegan","rafalib","parallel")
inst_packages(pkgs)
#---------



#############################
### LOAD Seurat.v3 OBJECT ###
#############################
cat("\nLoading/ data and metadata ...\n")
DATA <- readRDS(opt$Seurat_object_path)
cat("The total dimensions of your dataset is: ",dim(DATA),"\n")
#---------



######################################################
### CALCULATE DIVERSITY INDEXES OF GENE EXPRESSION ###
######################################################
cat("\nCalculating data diveristy indexes ...\n")
indexes <- t(apply(DATA@assays[[opt$assay]]@counts,2,function(x) {
      c(vegan::diversity(x,index = "simpson"),
       vegan::diversity(x,index = "invsimpson"),
       vegan::diversity(x,index = "shannon"),
       Gini(x)) }))
DATA$simp_index <- indexes[,1]
DATA$invsimp_index <- indexes[,2]
DATA$shan_index <- indexes[,3]
DATA$gini_index <- indexes[,4]
invisible(gc())
#---------



#############################################
### CALCULATE PERCENTAGE OF GENE FAMILIES ###
#############################################
cat("\nCalculating percentage of mitocondrial/ribosomal genes ...\n")
Gene.groups <- substring(rownames(x = DATA@assays[[opt$assay]]@counts),1,3)
seq_depth <- Matrix::colSums(DATA@assays[[opt$assay]]@counts)
temp <- rowsum(as.matrix(DATA@assays[[opt$assay]]@counts),Gene.groups)
perc <- sort(apply( t(temp) / seq_depth,2,median) ,decreasing = T)*100
tot <- sort(rowSums(temp)/sum(temp),decreasing = T)*100

#Compute the relative expression of each gene per cell
rel_expression <- Matrix::t( Matrix::t(DATA@assays[[opt$assay]]@counts) / Matrix::colSums(DATA@assays[[opt$assay]]@counts)) * 100
most_expressed <- sort(apply(rel_expression,1,median),T)[1:100] / ncol(DATA)

png(filename = paste0(opt$output_path,"/Gene_familty proportions.png"),width = 600*3,height = 4*600,res = 150)
mypar(4,1,mar=c(5,5,2,1))
boxplot( as.matrix(Matrix::t(rel_expression[names(most_expressed),])),cex=.1,outline=T,las=2,main="% total count per cell",col=hue_pal()(100))
boxplot( (t(temp)/seq_depth) [,names(perc)[1:100]]*100,outline=T,las=2,main="% reads per cell",col=hue_pal()(100))
boxplot(t(temp)[,names(perc)[1:100]], outline=T,las=2,main="reads per cell",col=hue_pal()(100) )
barplot(tot[names(tot)[1:100]],las=2,xaxs="i",main="Total % reads (all cells)",col=hue_pal()(100))
invisible(dev.off())

for(i in unique( c("rpl","rps","mito",unlist(strsplit(casefold(opt$plot_gene_family),","))))){
  cat(i,"\t")
  family.genes <- rownames(DATA@assays[[opt$assay]]@counts)[grep(pattern = paste0("^",ifelse(i=="mito","mt-",i)), x = casefold(rownames(DATA@assays[[opt$assay]]@counts)), value = F)]
  if(length(family.genes)>1){DATA <- PercentageFeatureSet(DATA,features = family.genes,assay = opt$assay,col.name = paste0("perc_",ifelse(i=="mt-","mito",i)) )}
}

rm("temp","perc","tot","Gene.groups","i","indexes")
invisible(gc())
#---------



############################################
### CALCULATING GENE BIOTYPE PERCENTAGES ###
############################################
cat("\nCalculating gene biotype percentages ...\n")
mart = useMart("ensembl", dataset = paste0(opt$species_use,"_gene_ensembl"),host="apr2019.archive.ensembl.org")
annot <- getBM(c("external_gene_name","gene_biotype"),mart = mart)
gene_biotype <- annot[match(rownames(DATA@assays[[opt$assay]]@counts) , annot[,1]),2]
gene_biotype[is.na(gene_biotype)] <- "unknown"

png(filename = paste0(opt$output_path,"/Gene_biotype_proportions.png"),width = 600*3,height = 600,res = 150)
mypar(1,3,mar=c(4,2,2,1))
pie(sort(table(gene_biotype),decreasing = T), clockwise = T,col = hue_pal()(length(unique(gene_biotype))))
title("before filtering")
par(mar=c(10,2,2,1))
barplot(sort(table(gene_biotype),decreasing = T),las=2,xaxs="i",main="Total reads (all cells)",col=hue_pal()(100))

temp <- rowsum(as.matrix(DATA@assays[[opt$assay]]@counts),group=gene_biotype)
o <- order(apply(temp,1,median),decreasing = T)
boxplot( (t(temp)/Matrix::colSums(DATA@assays[[opt$assay]]@counts))[,o]*100,outline=F,las=2,main="% reads per cell",col=hue_pal()(100))
invisible(dev.off())

aaa <- setNames(as.data.frame(((t(temp)/Matrix::colSums(DATA@assays[[opt$assay]]@counts))[,o]*100)[,names(sort(table(gene_biotype),decreasing = T))[1:6]]),paste0("perc_",names(sort(table(gene_biotype),decreasing = T))[1:6]))
DATA@meta.data <- cbind(DATA@meta.data,aaa)
#---------



########################
### NORMALIZING DATA ###
########################
cat("\nNormalizing counts ...\n")
#NOTE: Seurat.v3 has some issues with filtering, so we need to re-create the object for this step
DATA <- NormalizeData(object = DATA, scale.factor = 1000)
#---------



##########################
### CELL CYCLE SCORING ###
##########################
cat("\nPredicting cell cycle scores with Seurat ...\n")
s.genes <- Seurat::cc.genes$s.genes
g2m.genes <- Seurat::cc.genes$g2m.genes

if(casefold(opt$species_use) != "hsapiens"){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="apr2019.archive.ensembl.org")
  mart = useMart("ensembl", dataset = paste0(casefold(opt$species_use),"_gene_ensembl"),host="apr2019.archive.ensembl.org" )
  s.genes = getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = s.genes , mart = mart, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=F,valuesL = "hgnc_symbol")[,1]
  g2m.genes = getLDS(attributes = c("external_gene_name"), filters = "external_gene_name", values = g2m.genes , mart = mart, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=F,valuesL = "hgnc_symbol")[,1]
}

DATA <- CellCycleScoring(object = DATA, s.features = s.genes, g2m.features = g2m.genes)
DATA$G1.Score <- 1 - ( DATA$S.Score + DATA$G2M.Score )
DATA$CC.Diff <- DATA$S.Score - DATA$G2M.Score
#---------



###############
### PLOT QC ###
###############
cat("\nPlotting QC metrics ...\n")
for(i in as.character(unlist(strsplit(opt$columns_metadata,",")))){
feats <- colnames(DATA@meta.data) [ grepl("nFeature|nCount|perc|_index|[.]Score",colnames(DATA@meta.data) ) ]
png(filename = paste0(opt$output_path,"/QC_",i,"_ALL.png"),width = 1200*(length(unique(DATA@meta.data[,i]))/2+1),height = 700*ceiling(length(feats)/5),res = 200)
print(VlnPlot(object = DATA, features  = feats, ncol = 5,group.by = i,pt.size = .1,assay = opt$assay))
invisible(dev.off())}
#---------




#######################################
### SAVING THE RAW Seurat.v3 OBJECT ###
#######################################
cat("\nNumber of cells per metadata parameter for RAW Unfiltered data...\n")
for(i in strsplit(opt$columns_metadata,",")[[1]] ){
  cat("\n",i)
  print(table( DATA@meta.data[,i] ))
}
cat("\nDimentions of the raw.data objects BEFORE filtering ...\n")
print( dim(DATA@assays[[opt$assay]]@counts) )

cat("\nSaving the RAW Seurat object ...\n")
write.csv2(DATA@meta.data,paste0(opt$output_path,"/QC_metadata_all_cells.csv"),row.names = T)
saveRDS(DATA, file = paste0(opt$output_path,"/raw_seurat_object.rds") )
#---------



###########################################
### SELECTING ONLY PROTEIN-CODING GENES ###
###########################################
cat("\nSelect only the protein-coding genes ...\n")
if( casefold(opt$remove_non_coding) == 'true' ){
  sel <- annot[match(rownames(DATA@assays[[opt$assay]]@counts) , annot[,1]),2] == "protein_coding"
  genes_use <- rownames(DATA@assays[[opt$assay]]@counts)[sel]
  genes_use <- as.character(na.omit(genes_use))
  DATA@assays[[opt$assay]]@counts <- DATA@assays[[opt$assay]]@counts[genes_use,]
}
#---------




#############################################
### REMOVING SELECTED GENES FROM THE DATA ###
#############################################
cat("\nRemoving selected genes from the data ...\n")
print( strsplit(opt$remove_gene_family,",")[[1]] )
if(opt$remove_gene_family != "none"){
  genes_use <- rownames(DATA@assays[[opt$assay]]@counts)[!grepl(gsub(",","|",casefold(opt$remove_gene_family) ) , casefold(rownames(DATA@assays[[opt$assay]]@counts)))]
  DATA@assays[[opt$assay]]@counts <- DATA@assays[[opt$assay]]@counts[genes_use,]
}
#---------



######################
### CELL FILTERING ###
######################
cat("\nFiltering low quality cells ...\n")
NF <-  DATA@meta.data [ grepl("nFeature",colnames(DATA@meta.data)) ][,1]
NC <-  DATA@meta.data [ grepl("nCount",colnames(DATA@meta.data)) ][,1]

Ts <- data.frame(
  MitoT = between(DATA$perc_mito,0.00,25),
  RpsT = between(DATA$perc_rps,3,50),
  RplT = between(DATA$perc_rps,3,50),
  nUMIT = between(NF,quantile(NF,probs = c(0.005)),quantile(NF,probs = c(0.995))),
  nCountT = between(NC,quantile(NC,probs = c(0.005)),quantile(NC,probs = c(0.995))),
  GiniT = between(DATA$gini_index,0.9,1),
  SimpT = between(DATA$simp_index,0.95,1),
  protein_codingT = between(DATA$perc_protein_coding,0.8,1),
  row.names = rownames(DATA@meta.data) )
print(head(Ts,20))

dim(DATA)
cell_use <- rownames(Ts)[ rowSums(!Ts) == 0 ]
length(cell_use)
#---------



####################################
### RE-NORMALIZING FILTERED DATA ###
####################################
cat("\nNormalizing counts ...\n")
#NOTE: Seurat.v3 has some issues with filtering, so we need to re-create the object for this step
DATA <- CreateSeuratObject(counts = DATA@assays[[opt$assay]]@counts[,cell_use] , assay = opt$assay, meta.data = DATA@meta.data[cell_use,], min.cells = as.numeric(opt$min_gene_count),min.features = as.numeric(opt$min_gene_per_cell))
DATA <- NormalizeData(object = DATA,scale.factor = 1000)

cat("\nDimentions of the raw.data objects AFTER filtering ...\n")
print( dim(DATA@assays[[opt$assay]]@counts) )
#---------



###############
### PLOT QC ###
###############
cat("\nPlotting QC metrics ...\n")
for(i in as.character(unlist(strsplit(opt$columns_metadata,",")))){
png(filename = paste0(opt$output_path,"/QC_",i,"_FILTERED.png"),width = 1200*(length(unique(DATA@meta.data[,i]))/2+1),height = 700*ceiling(length(feats)/5),res = 200)
print(VlnPlot(object = DATA, features  = feats, ncol = 5,group.by = i,pt.size = .1,assay = opt$assay))
invisible(dev.off())}
#---------



########################################
### QUANTIFICATION OF DETECTED GENES ###
########################################
# cat("\nIdentification of detected genes ...\n")
# library(pheatmap)
# rawdata <- as.matrix(DATA@assays[[opt$assay]]@counts)[!(rowSums(as.matrix(DATA@assays[[opt$assay]]@counts)) == 0),]
# N <- ncol(rawdata)
# filter_test <- data.frame("Exp>1 in 5 cells"=(rowSums(rawdata >= 1) >= 5)*1,
#                           "Exp>1 in 10 cells"=(rowSums(rawdata >= 1) >= 10)*1,
#                           "Exp>1 in 20 cells"=(rowSums(rawdata >= 1) >= 20)*1,
#                           "Exp>1 in 30 cells"=(rowSums(rawdata >= 1) >= 30)*1,
#                           "Exp>1 in 40 cells"=(rowSums(rawdata >= 1) >= 40)*1,
#                           "Exp>1 in 60 cells"=(rowSums(rawdata >= 1) >= 60)*1,
#                           "Exp>3 in 2 cells"=(rowSums(rawdata >= 3) >= 2)*1,
#                           "Exp>3 in 5 cells"=(rowSums(rawdata >= 3) >= 5)*1,
#                           "Exp>3 in 10 cells"=(rowSums(rawdata >= 3) >= 10)*1,
#                           "Exp>3 in 20 cells"=(rowSums(rawdata >= 3) >= 20)*1,
#                           "Exp>5 in 1 cell"=(rowSums(rawdata >= 5) >= 1)*1,
#                           "Exp>5 in 2 cells"=(rowSums(rawdata >= 5) >= 2)*1,
#                           "Exp>5 in 3 cells"=(rowSums(rawdata >= 5) >= 3)*1,
#                           "Exp>5 in 5 cells"=(rowSums(rawdata >= 5) >= 5)*1,
#                           "RowMeans >0.01"=(rowMeans(rawdata) >= 0.01)*1,
#                           "RowMeans >0.02"=(rowMeans(rawdata) >= 0.02)*1,
#                           "RowMeans >0.05"=(rowMeans(rawdata) >= 0.05)*1,
#                           "RowMeans >0.10"=(rowMeans(rawdata) >= 0.10)*1,
#                           "RowMeans >0.30"=(rowMeans(rawdata) >= 0.30)*1,
#                           "Exp>1 in 1pct cells"=(rowSums(rawdata >= 1) >= N/100)*1,
#                           "Exp>1 in 2.5pct cells"=(rowSums(rawdata >= 1) >= N/100*2.5)*1,
#                           "Exp>1 in 5pct cells"=(rowSums(rawdata >= 1) >= N/100*5)*1,
#                           "Exp>1 in 7.5pct cells"=(rowSums(rawdata >= 1) >= N/100*7.5)*1,
#                           "Exp>1 in 10pct cells"=(rowSums(rawdata >= 1) >= N/10)*1)
# pheatmap(t(filter_test),color=c("grey80","navy"),cluster_rows = F,gaps_row = c(6,10,14,19),
#          filename = paste0(opt$output_path,"/Heatmap_min_gene_count.png"))
# 
# png(filename = paste0(opt$output_path,"/Barplot_min_gene_count.png"),width = 1500,height = 1000,res = 200)
# par(mar=c(10,4,2,1))
# barplot(colSums(filter_test),las=2,yaxs="i",border=NA)
# invisible(dev.off())
# rm("rawdata")
#---------



#####################################
### SAVING FILTERED SEURAT OBJECT ###
#####################################
cat("\nNumber of cells per metadata parameter for raw FILTERED data...\n")
for(i in strsplit(opt$columns_metadata,",")[[1]] ){
  cat("\n",i)  ;   print(table( DATA@meta.data[,i] )) }

cat("\nSaving filtered Seurat object ...\n")
saveRDS(DATA, file = paste0(opt$output_path,"/filt_seurat_object.rds") )
#---------



#############################
### SYSTEM & SESSION INFO ###
#############################
#---------
print_session_info()
#---------
