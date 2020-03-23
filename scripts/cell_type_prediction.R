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
  make_option(c("-c", "--clustering_use"),        type = "character",   metavar="character",   default='none',   help="The clustering name to be used for the analysis"),
  make_option(c("-a", "--assay"),                 type = "character",   metavar="character",   default='RNA',  help="Assay to be used in the analysis."),
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
suppressMessages(suppressWarnings({
  library(Seurat)
  library(scales)
  library(fields)
  library(data.table)
  library(ggplot2)
}))

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
marker_lists <- unlist(strsplit(opt$marker_lists,","))
print(marker_lists)

for(i in sub(".*/","",sub(".csv","",marker_lists)) ){
  PATH <- paste0(opt$output_path,"/",sub(".*/","",i))
  if(!dir.exists(PATH)){dir.create(PATH,recursive = T)}
  
  cat("\nProcessing list '", sub(".*/","",i) ,"' ...\n")
  cellIDs <- read.csv2(marker_lists[grep(i,marker_lists)],header =T)
  cellIDs <- as.list(as.data.frame(cellIDs))
  cellIDs <- lapply(cellIDs, function(x) casefold( as.character(x[x!=""]) ) )
  #cellIDs <- lapply(cellIDs, function(x) x[1:min(10,length(x))] )
  print(cellIDs)
  

cat("\nCreating cell identity matrix\n")
cell_ident <- unique(unlist(cellIDs))
#cell_ident <- lapply(cellIDs,function(x) cell_ident %in% x)
cell_ident <- lapply(cellIDs,function(x) ifelse(cell_ident %in% x,1,ifelse(paste0("-",cell_ident) %in% x,-1,0)))
cell_ident <- as.data.frame(cell_ident)
rownames(cell_ident) <- casefold(unique(unlist(cellIDs)))
#print(cell_ident[1:10,1:10])

#hist(rowSums(cell_ident))
#--------------------------------



#Filter data for the genes inthe cell identity matrix
#--------------------------------
cat("\nSelecting detected genes ...\n")
sel <- rownames(DATA@assays[[opt$assay]]@data) [ casefold(rownames(DATA@assays[[opt$assay]]@data)) %in% rownames(cell_ident) ]
#sel <- rownames(cell_ident)[ rownames(cell_ident) %in% casefold(rownames(DATA@assays[[opt$assay]]@data)) ]
cell_ident <- cell_ident[ casefold(sel ) ,]
gc()
#--------------------------------



#Using correlations as distance
#--------------------------------
cat("\nPredicting cell type by correlation method ...\n")

cat("\nComputing correlations ...\n")
cors <- apply(DATA@assays[[opt$assay]]@data[ sel ,],2,function(x) cor(x , cell_ident) )
cors[is.na(cors)] <- -1
rownames(cors) <- colnames(cell_ident)
cors2 <- t(t(cors) / apply(cors,2,max))
cors2[1:nrow(cors2),1:20]
print(cors2[,1:5])
gc()
try(write.csv2(cors,paste0(opt$output_path,"/",i,"/cell_pred_correlation_",i,".csv"),row.names = T))

cat("\nPredicting cell types ...\n")
pred <- unlist( apply(cors2,2,function(x) colnames(cell_ident) [which.max(x)]) )
my_nas <- colnames(cors2)[! colnames(cors2) %in% names(pred)]
pred <- c(pred , setNames(rep(NA,length(my_nas)),my_nas))

cat("\nPlotting ...\n")
DATA <- AddMetaData(DATA,metadata = pred,col.name = paste0("cell_pred_correlation_",i))

for(j in names(DATA@reductions) ){
  temp2 <- DimPlot(DATA,dims = 1:2,reduction = j,group.by = paste0("cell_pred_correlation_",i), pt.size = .3,ncol = 8)
  ggplot2::ggsave(temp2,filename = paste0("cell_cluster_pred_correlation_",j,".png"), path = paste0(opt$output_path,"/",i), dpi = 300,units = "mm",width = 140,height = 110,limitsize = FALSE )
}

#png(filename = paste0(opt$output_path,"/",i,"/tSNE_cell_cluster_pred_correlation.png"),width = 700,height = 600,res = 150)
#TSNEPlot(object = DATA,group.by=paste0("cell_pred_correlation_",i),pt.size = .3,plot.title= paste0("cell_pred_correlation_",i))
#invisible(dev.off())

png(filename = paste0(opt$output_path,"/",i,"/UMAP_single_cell_pred_correlation.png"),width = 400*4,height = 350*ceiling(nrow(cors) / 4),res = 150)
par(mar=c(1.5,1.5,3,5), mfrow=c(ceiling(nrow(cors) / 4),4))
myCorGrad <- paste0(colorRampPalette(c("gray85","gray85","gray70","orange3","firebrick","red"))(10))
for(j in rownames(cors)){
  lim <- max(as.numeric(cors[j,]),na.rm = T)
  temp <- ((cors[j,]) - 0) / ( max(lim,0.6) + 0)
  temp[is.na(temp)] <- min(temp,na.rm = T)
  temp <- round((temp)*9)+1
  temp[temp <= 1] <- 1
  o <- order(temp)
  plot(DATA@reductions[["umap"]]@cell.embeddings[o,],pch=20,cex=0.6, line=0.5, col=myCorGrad[ temp[o] ], yaxt="n",xaxt="n",xlab="tSNE1",ylab="tSNE2",lwd=0.25, main=paste0("Cor. to ",j))
  image.plot(1,1,cbind(0,lim),legend.only = T,col = myCorGrad)
}
dev.off()
#--------------------------------




#Using normalized euclidean distance
#--------------------------------
cat("\nPredicted cell type by euclidean-distance method ...\n")
c <- apply(DATA@assays[[opt$assay]]@data[ sel ,],2,function(x) cell_ident - (x/max(x)) )
c1 <- lapply(c,function(x) colSums(x^2) )
c2 <- as.data.frame(c1)
print(c2[1:10,1:10])

write.csv(c2,paste0(opt$output_path,"/",i,"/cell_pred_euclidean_",i,".csv"),row.names = T)


pred2 <- unlist(apply(c2,2,function(x) colnames(cell_ident) [ which.min(x) ]))
my_nas <- colnames(cors)[! colnames(cors) %in% names(pred2)]
#pred2 <- c(pred2 , setNames(rep(NA,length(my_nas)),my_nas))

cat("\nPlotting ...\n")
#DATA@meta.data[[paste0("cell_pred_euclidean_",i)]] <- c(pred2)
DATA <- AddMetaData(DATA,metadata = factor(pred2), col.name = paste0("cell_pred_euclidean_",i))
for(j in names(DATA@reductions) ){
  temp2 <- DimPlot(DATA,dims = 1:2,reduction = j,group.by = paste0("cell_pred_euclidean_",i), pt.size = .3,ncol = 8)
  ggplot2::ggsave(temp2,filename = paste0("cell_cluster_pred_euclidean_",j,".png"), 
                  path = paste0(opt$output_path,"/",i), dpi = 300,units = "mm",width = 140,height = 110,limitsize = FALSE )
}

#png(filename = paste0(opt$output_path,"/",i,"/tSNE_cell_cluster_pred_euclidean.png"),width = 700,height = 600,res = 150)
#TSNEPlot(object = DATA,group.by=paste0("cell_pred_euclidean_",i),pt.size = .3,plot.title= paste0("cell_pred_euclidean_",i))
#invisible(dev.off())

png(filename = paste0(opt$output_path,"/",i,"/UMAP_single_cell_pred_euclidean.png"),width = 400*4,height = 350*ceiling(nrow(cors) / 4),res = 150)
par(mar=c(1.5,1.5,3,5), mfrow=c(ceiling(nrow(c2) / 4),4))
for(j in rownames(c2)){
  myCorGrad <- paste0(colorRampPalette(c("red","firebrick","orange3","gray70","gray70","gray85","gray85","gray85"))(10))
  lim <- max(as.numeric(c2[j,]),na.rm = T)
  temp <- ((as.numeric(c2[j,]) ) + 0) / (lim + 0)
  temp[is.na(temp)] <- min(temp,na.rm = T)
  temp <- round((temp)*9)+1
  temp[temp <= 1] <- 1
  o <- order(temp,decreasing = T)
  plot(DATA@reductions[["umap"]]@cell.embeddings[o,],pch=20,cex=0.6, line=0.5, col=myCorGrad[ temp[o] ], yaxt="n",xaxt="n",xlab="tSNE1",ylab="tSNE2",lwd=0.25, main=paste0("Cor. to ",j))
  image.plot(1,1,cbind(0,lim),legend.only = T,col = myCorGrad)
}
dev.off()

#--------------------------------
}


#########################################################
### Plot the percentage of cell types in each cluster ###
#########################################################
if (!(casefold(opt$clustering_use) == 'none')){
  cat("\nPlotting cell type distribution across clusters...\n")
  
  DATA <- SetIdent(DATA, value = factor(as.character(DATA@meta.data[,opt$clustering_use])))
  
  pred_factor <- as.factor(DATA@meta.data[, paste0("cell_pred_correlation_",i)])
  cell_type_levels <- levels(pred_factor)
  
  proportion <- as.data.frame(lapply(levels(DATA@active.ident),function(x){c(unname(table(pred_factor[DATA@active.ident==x]))) [1:length(cell_type_levels)]} ))
  proportion[is.na(proportion)] <- 0
  rownames(proportion) <- cell_type_levels
  colnames(proportion) <- levels(DATA@active.ident)
  
  cl_order <- order(as.numeric(levels(DATA@active.ident)))
  sa <- cbind(stack(as.data.frame(t(t(proportion[,cl_order])/colSums(proportion[,cl_order])) )), rep(rownames(proportion),ncol(proportion)) )
  colnames(sa) <- c("prop","clusters","celltype")
  
  png(filename = paste0(opt$output_path,"/",i,"/cell_cluster_pred_correlation_barplot.png"),width = 1100,height = 500,res = 150)
  print(ggplot(data=sa, aes(x=clusters, y=prop, fill=celltype) ) + geom_col())
  invisible(dev.off())
}
#--------- 



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



