#!/usr/bin/env Rscript

### LOAD OPTPARSE
#---------
if(!require("optparse")){install.packages("optparse",, repos='http://cran.us.r-project.org')};
library(optparse)
#---------

### DEFINE PATH TO LOCAL FILES
#---------
cat("\nRunning QUALITY CONTROL with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object"),
  make_option(c("-c", "--columns_metadata"),      type = "character",   metavar="character",   default='none',  help="Column names in the Metadata matrix (only factors allowed, not continuous variables)"),
  make_option(c("-s", "--species_use"),           type = "character",   metavar="character",   default='none',  help="Species from the sample for cell scoring"),
  make_option(c("-p", "--cell_phase_info"),       type = "character",   metavar="character",   default='none',  help="Path for the cell cycle phase genes"),
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



### LOAD Seurat OBJECT 
#---------
cat("\nLoading/ data and metadata ...\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------




### Calculate the distribution of detected genes in each sample
#---------
cat("\nCalculating data diveristy indexes ...\n")
gini_ind <- apply(DATA@raw.data,2,Gini)
DATA <- AddMetaData(object = DATA, metadata = gini_ind, col.name = "gini_ind")

simp_ind <- apply(DATA@raw.data,2,function(x) vegan::diversity(x,index = "simpson"))
DATA <- AddMetaData(object = DATA, metadata = simp_ind, col.name = "simp_ind")

invsimp_ind <- apply(DATA@raw.data,2,function(x) vegan::diversity(x,index = "invsimpson"))
DATA <- AddMetaData(object = DATA, metadata = invsimp_ind, col.name = "invsimp_ind")

shan_ind <- apply(DATA@raw.data,2,function(x) vegan::diversity(x,index = "shannon"))
DATA <- AddMetaData(object = DATA, metadata = shan_ind, col.name = "shan_ind")
#---------



### Calculate the percentage of counts from structural genes (mitocondrial and ribosomal)
#---------
cat("\nCalculating percentage of mitocondrial/ribosomal genes ...\n")
Gene.groups <- substring(rownames(x = DATA@data),1,3)
temp <- rowsum(as.matrix(DATA@raw.data),Gene.groups) / colSums(as.matrix(DATA@raw.data))
perc <- sort(rowMeans(temp),decreasing = T)

png(filename = paste0(opt$output_path,"/Gene_familty proportions.png"),width = 400*3,height = 2*400,res = 150)
mypar(2,1)
boxplot(100*t(temp[rownames(temp)%in%names(perc)[1:40],])[,names(perc)[1:40]],outline=F,las=2,ylab="% reads per cell",col=hue_pal()(40) )
barplot(perc[1:40]*100,las=2,xaxs="i",ylab="mean % reads",col=hue_pal()(40))
dev.off()


mito.genes <- grep(pattern = "^mt-", x = rownames(x = DATA@data), value = TRUE,ignore.case = T)
percent.mito <- Matrix::colSums(DATA@raw.data[mito.genes, ]) / Matrix::colSums(DATA@raw.data)
DATA <- AddMetaData(object = DATA, metadata = percent.mito, col.name = "percent.mito")

Rps.genes <- grep(pattern = "^Rps[123456789]", x = rownames(x = DATA@data), value = TRUE,ignore.case = T)
percent.Rps <- Matrix::colSums(DATA@raw.data[Rps.genes, ]) / Matrix::colSums(DATA@raw.data)
DATA <- AddMetaData(object = DATA, metadata = percent.Rps, col.name = "percent.Rps")

Rpl.genes <- grep(pattern = "^Rpl[123456789]", x = rownames(x = DATA@data), value = TRUE,ignore.case = T)
percent.Rpl <- Matrix::colSums(DATA@raw.data[Rpl.genes, ]) / Matrix::colSums(DATA@raw.data)
DATA <- AddMetaData(object = DATA, metadata = percent.Rpl, col.name = "percent.Rpl")
#---------



### Plotting QC plots
#---------
for(i in as.character(unlist(strsplit(opt$columns_metadata,",")))){
  png(filename = paste0(opt$output_path,"/QC_",i,".png"),width = 1200*(length(unique(DATA@meta.data[,i]))/2+1),height = 1400,res = 200)
  print(VlnPlot(object = DATA, features.plot = c("nGene", "nUMI", "percent.mito","percent.Rps","percent.Rpl","shan_ind","simp_ind","gini_ind","invsimp_ind"), nCol = 5,group.by = i,point.size.use = .1))
  dev.off()}
#---------




### Normalizing the data
#---------
cat("\nNormalizing counts ...\n")
DATA <- NormalizeData(object = DATA)
#---------


### Identification of cell cycle phase using Seurat
#---------
cat("\nPredicting cell cycle scores with Seurat ...\n")
cc.genes <- readLines(con = paste0(opt$cell_phase_info,"/regev_lab_cell_cycle_genes.txt"))
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

if(opt$species_use == "mouse"){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  s.genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = s.genes , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=F,valuesL = "hgnc_symbol")[,1]
  g2m.genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = g2m.genes , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=F,valuesL = "hgnc_symbol")[,1]
}

DATA <- CellCycleScoring(object = DATA, s.genes = s.genes, g2m.genes = g2m.genes)
DATA@meta.data$CC.Diff <- DATA@meta.data$S.Score - DATA@meta.data$G2M.Score
#---------





### Saving the RAW Seurat object
#---------
cat("\nSaving the RAW Seurat object ...\n")
write.csv(DATA@meta.data,paste0(opt$output_path,"/QC_metadata_all_cells.csv"),row.names = T)
saveRDS(DATA, file = paste0(opt$output_path,"/Raw_Seurat_Object.rds") )
#---------



### Cell filtering
#---------
cat("\nFiltering low quality cells ...\n")
Ts <- data.frame(
nGeneT = between(DATA@meta.data$nGene,quantile(DATA@meta.data$nGene,probs = c(.01)),quantile(DATA@meta.data$nGene,probs = c(.99))),
MitoT = between(DATA@meta.data$percent.mito,quantile(DATA@meta.data$percent.mito,probs = c(.01)),quantile(DATA@meta.data$percent.mito,probs = c(.99))),
nUMIT = between(DATA@meta.data$nUMI,quantile(DATA@meta.data$nUMI,probs = c(.01)),quantile(DATA@meta.data$nUMI,probs = c(.99))),
RpsT = between(DATA@meta.data$percent.Rps,quantile(DATA@meta.data$percent.Rps,probs = c(.01)),quantile(DATA@meta.data$percent.Rps,probs = c(.99))),
RplT = between(DATA@meta.data$percent.Rpl,quantile(DATA@meta.data$percent.Rpl,probs = c(.01)),quantile(DATA@meta.data$percent.Rpl,probs = c(.99))),
ShanT = between(DATA@meta.data$shan_ind,quantile(DATA@meta.data$shan_ind,probs = c(.01)),quantile(DATA@meta.data$shan_ind,probs = c(.99))),
row.names = rownames(DATA@meta.data)
)

DATA <- SubsetData(DATA,cells.use = rownames(Ts)[rowSums(!Ts) == 0])

for(i in as.character(unlist(strsplit(opt$columns_metadata,",")))){
  png(filename = paste0(opt$output_path,"/QC_",i,"_FILTERED.png"),width = 1200*(length(unique(DATA@meta.data[,i]))/2+1),height = 1400,res = 200)
  print(VlnPlot(object = DATA, features.plot = c("nGene", "nUMI", "percent.mito","percent.Rps","percent.Rpl","shan_ind","simp_ind","gini_ind","invsimp_ind"), nCol = 5,group.by = i,point.size.use = .1))
  dev.off()}
#---------




### Identification of detected genes
#---------
cat("\nIdentification of detected genes ...\n")
library(pheatmap)
rawdata <- as.matrix(DATA@raw.data)[!(rowSums(as.matrix(DATA@raw.data)) == 0),]
N <- ncol(rawdata)
filter_test <- data.frame("Exp>1 in 5 cells"=(rowSums(rawdata >= 1) >= 5)*1,
                          "Exp>1 in 10 cells"=(rowSums(rawdata >= 1) >= 10)*1,
                          "Exp>1 in 20 cells"=(rowSums(rawdata >= 1) >= 20)*1,
                          "Exp>1 in 30 cells"=(rowSums(rawdata >= 1) >= 30)*1,
                          "Exp>1 in 40 cells"=(rowSums(rawdata >= 1) >= 40)*1,
                          "Exp>1 in 60 cells"=(rowSums(rawdata >= 1) >= 60)*1,
                          "Exp>3 in 2 cells"=(rowSums(rawdata >= 3) >= 2)*1,
                          "Exp>3 in 5 cells"=(rowSums(rawdata >= 3) >= 5)*1,
                          "Exp>3 in 10 cells"=(rowSums(rawdata >= 3) >= 10)*1,
                          "Exp>3 in 20 cells"=(rowSums(rawdata >= 3) >= 20)*1,
                          "Exp>5 in 1 cell"=(rowSums(rawdata >= 5) >= 1)*1,
                          "Exp>5 in 2 cells"=(rowSums(rawdata >= 5) >= 2)*1,
                          "Exp>5 in 3 cells"=(rowSums(rawdata >= 5) >= 3)*1,
                          "Exp>5 in 5 cells"=(rowSums(rawdata >= 5) >= 5)*1,
                          "RowMeans >0.01"=(rowMeans(rawdata) >= 0.01)*1,
                          "RowMeans >0.02"=(rowMeans(rawdata) >= 0.02)*1,
                          "RowMeans >0.05"=(rowMeans(rawdata) >= 0.05)*1,
                          "RowMeans >0.10"=(rowMeans(rawdata) >= 0.10)*1,
                          "RowMeans >0.30"=(rowMeans(rawdata) >= 0.30)*1,
                          "Exp>1 in 1pct cells"=(rowSums(rawdata >= 1) >= N/100)*1,
                          "Exp>1 in 2.5pct cells"=(rowSums(rawdata >= 1) >= N/100*2.5)*1,
                          "Exp>1 in 5pct cells"=(rowSums(rawdata >= 1) >= N/100*5)*1,
                          "Exp>1 in 7.5pct cells"=(rowSums(rawdata >= 1) >= N/100*7.5)*1,
                          "Exp>1 in 10pct cells"=(rowSums(rawdata >= 1) >= N/10)*1)
pheatmap(t(filter_test),color=colorRampPalette(c("grey80","grey80","navy"))(99),filename = paste0(opt$output_path,"/Heatmap_gene_filtering.png"))

png(filename = paste0(opt$output_path,"/Barplot_gene_filtering.png"),width = 1500,height = 1000,res = 200)
par(mar=c(10,4,2,1))
barplot(sort(colSums(filter_test),decreasing = T),las=2,yaxs="i",border=NA)
dev.off()
#---------







### Saving the Seurat object
#---------
cat("\nSaving filtered Seurat object ...\n")
saveRDS(DATA, file = paste0(opt$output_path,"/Filt_Seurat_Object.rds") )
#---------


cat("\n!!! Script executed Sucessfully !!!\n")



### System and session information
#---------
cat("\n\n\n\n... SYSTEM INFORMATION ...\n")
Sys.info()

cat("\n\n\n\n... SESSION INFORMATION ...\n")
sessionInfo()
#---------