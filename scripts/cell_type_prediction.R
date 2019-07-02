#!/usr/bin/env Rscript

### LOAD OPTPARSE
#---------
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')};
library(optparse)
#---------


### DEFINE PATH TO LOCAL FILES
#---------
cat("\nRunning CELL TYPE PREDICTION with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object"),
  make_option(c("-m", "--marker_lists"),          type = "character",   metavar="character",   default='none',  help="A folder containing .csv files with the list of markers for comparison"),
  make_option(c("-c", "--cluster_use"),           type = "character",   metavar="character",   default='all',   help="The clustering name to be used for the analysis"),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output directory")
) 
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)
#---------



### DEFINE PATH TO LOCAL FILES
#---------
col_scale <- c("grey85","navy")
myCorGrad <- colorRampPalette(c("gray85","gray85","gray70","orange3","firebrick","red"))(9)
#if(!dir.exists(paste0(opt$output_path,"/cell_prediction"))){dir.create(paste0(opt$output_path,"/cell_prediction"),recursive = T)}
#---------



### LOAD LIBRARIES
#---------
cat("\nLoading/installing libraries ...\n")
initial.options <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(sub("--file=","",initial.options[grep("--file=",initial.options)]))
source( paste0(script_path,"/inst_packages.R") )
pkgs <- c("Seurat","scales","fields","data.table")
inst_packages(pkgs)
#---------




### LOAD Seurat OBJECT 
#---------
cat("\nLoading/ data and metadata ...\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------






### LOAD cell marker list
#---------
cat("\nLoading cell marker lists and creating cell identity matrix ...\n")
cat("\nThe following datasets were found ...\n")
marker_lists <- list.files( opt$marker_lists , pattern = ".csv")
print(marker_lists)

for(i in marker_lists){
  PATH <- paste0(opt$output_path,"/",i)
  if(!dir.exists(PATH)){dir.create(PATH,recursive = T)}
  

cat("\nProcessing list ", i ," ...\n")
cellIDs <- read.csv2(paste0(opt$marker_lists,"/",i),header = F)
rownames(cellIDs) <- as.character(cellIDs[,1])
print(cellIDs)
cellIDs <- as.list(as.data.frame(t(cellIDs[,-1])))
cellIDs <- lapply(cellIDs, function(x) as.character(x[x!=""]) )




#Create cell identity matrix
cell_ident <- unique(unlist(cellIDs))
#cell_ident <- lapply(cellIDs,function(x) cell_ident %in% x)
cell_ident <- lapply(cellIDs,function(x) ifelse(cell_ident %in% x,1,ifelse(paste0("-",cell_ident) %in% x,-1,0)))
cell_ident <- as.data.frame(cell_ident)
rownames(cell_ident) <- unique(unlist(cellIDs))
#--------------------------------



#Filter data for the genes inthe cell identity matrix
#--------------------------------
data <- as.matrix(DATA@data)
sel <- rownames(cell_ident)[ rownames(cell_ident) %in% rownames(data) ]
data <- data[sel,]
cell_ident <- cell_ident[sel,]
gc()
#--------------------------------




#Using correlations as distance
#--------------------------------
cat("\nPredicting cell type by correlation method ...\n")

cat("\nComputing correlations ...\n")
cors <- apply(data,2,function(x) cor(x , cell_ident) )
rownames(cors) <- colnames(cell_ident)
cors2 <- t(t(cors) / apply(cors,2,max))
cors2[1:nrow(cors2),1:20]

write.csv(cors,paste0(opt$output_path,"/Cell_pred_correlation_",i,".csv"),row.names = T)

cat("\nProdicting cell types ...\n")
pred <- unlist( apply(cors2,2,function(x) colnames(cell_ident) [which.max(x)]) )
my_nas <- colnames(cors2)[! colnames(cors2) %in% names(pred)]
pred <- c(pred , setNames(rep(NA,length(my_nas)),my_nas))

cat("\nPlotting ...\n")
DATA <- AddMetaData(DATA,metadata = pred,col.name = paste0("cell_pred_correlation_",i))
png(filename = paste0(opt$output_path,"/",i,"/tSNE_cell_cluster_pred_correlation.png"),width = 700,height = 600,res = 150)
TSNEPlot(object = DATA,group.by=paste0("cell_pred_correlation_",i),pt.size = .3,plot.title= paste0("cell_pred_correlation_",i))
invisible(dev.off())

png(filename = paste0(opt$output_path,"/",i,"/tSNE_single_cell_pred_correlation.png"),width = 400*4,height = 350*ceiling(nrow(cors) / 4),res = 150)
par(mar=c(1.5,1.5,3,5), mfrow=c(ceiling(nrow(cors) / 4),4))
myCorGrad <- paste0(colorRampPalette(c("gray85","gray85","gray70","orange3","firebrick","red"))(10))
for(j in rownames(cors)){
  lim <- max(as.numeric(cors[j,]),na.rm = T)
  temp <- ((cors[j,]) - 0) / ( max(lim,0.6) + 0)
  temp[is.na(temp)] <- min(temp,na.rm = T)
  temp <- round((temp)*9)+1
  temp[temp <= 1] <- 1
  o <- order(temp)
  plot(DATA@dr$tsne@cell.embeddings[o,],pch=20,cex=0.6, line=0.5, col=myCorGrad[ temp[o] ], yaxt="n",xaxt="n",xlab="tSNE1",ylab="tSNE2",lwd=0.25, main=paste0("Cor. to ",j))
  image.plot(1,1,cbind(0,lim),legend.only = T,col = myCorGrad)
}
dev.off()
#--------------------------------






#Using normalized euclidean distance
#--------------------------------
cat("\nPredicted cell type by euclidean-distance method ...\n")
c <- apply(data,2,function(x) cell_ident - (x/max(x)) )
c1 <- lapply(c,function(x) colSums(x^2) )
c2 <- as.data.frame(c1)
c2[1:nrow(c2),1:10]

write.csv(c2,paste0(opt$output_path,"/cell_pred_euclidean_",i,".csv"),row.names = T)


pred2 <- unlist(apply(c2,2,function(x) colnames(cell_ident) [x==min(x)]))
my_nas <- colnames(cors)[! colnames(cors) %in% names(pred2)]
pred2 <- c(pred2 , setNames(rep(NA,length(my_nas)),my_nas))

cat("\nPlotting ...\n")
DATA <- AddMetaData(DATA,metadata = pred2,col.name = paste0("cell_pred_euclidean_",i))
png(filename = paste0(opt$output_path,"/",i,"/tSNE_cell_cluster_pred_euclidean.png"),width = 700,height = 600,res = 150)
TSNEPlot(object = DATA,group.by=paste0("cell_pred_euclidean_",i),pt.size = .3,plot.title= paste0("cell_pred_euclidean_",i))
invisible(dev.off())

png(filename = paste0(opt$output_path,"/",i,"/tSNE_single_cell_pred_euclidean.png"),width = 400*4,height = 350*ceiling(nrow(cors) / 4),res = 150)
par(mar=c(1.5,1.5,3,5), mfrow=c(ceiling(nrow(c2) / 4),4))
for(j in rownames(c2)){
  myCorGrad <- paste0(colorRampPalette(c("red","firebrick","orange3","gray70","gray70","gray85","gray85","gray85"))(10))
  lim <- max(as.numeric(c2[j,]),na.rm = T)
  temp <- ((as.numeric(c2[j,]) ) + 0) / (lim + 0)
  temp[is.na(temp)] <- min(temp,na.rm = T)
  temp <- round((temp)*9)+1
  temp[temp <= 1] <- 1
  o <- order(temp,decreasing = T)
  plot(DATA@dr$tsne@cell.embeddings[o,],pch=20,cex=0.6, line=0.5, col=myCorGrad[ temp[o] ], yaxt="n",xaxt="n",xlab="tSNE1",ylab="tSNE2",lwd=0.25, main=paste0("Cor. to ",j))
  image.plot(1,1,cbind(0,lim),legend.only = T,col = myCorGrad)
}
dev.off()

#--------------------------------
}





### Saving the Seurat object
#---------
saveRDS(DATA, file = opt$Seurat_object_path )
#---------





cat("\n!!! Script executed Sucessfully !!!\n")


### System and session information
#---------
cat("\n\n\n\n... SYSTEM INFORMATION ...\n")
INFORMATION <- Sys.info()
print(as.data.frame(INFORMATION))

cat("\n\n\n\n... SESSION INFORMATION ...\n")
sessionInfo()
#---------



