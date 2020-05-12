#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
cat("\nRunning Gene module identification using WGCNA ... \n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object"),
  make_option(c("-c", "--columns_metadata"),      type = "character",   metavar="character",   default='none',  help="Column names in the Metadata matrix (only factors allowed, not continuous variables)"),
  make_option(c("-m", "--modules_of_interest"),   type = "character",   metavar="character",   default='none',  help="Gene modules of interest. After running this analysis, you can highlight some modules on the plot. Just provide the gene module numbers separated by a comma: '1,21,13,4,35,25'. Those will marked with an asterisk."),
  make_option(c("-g", "--gene_list_use"),         type = "character",   metavar="character",   default='none',  help="Path to a csv file containing which genes should be used for gene module identification. It is recommended to run differential expression and then use the list of differentially expressed genes as input here. A list of highlighly variabel genes is another alternative."),
  make_option(c("-k", "--number_of_modules"),     type = "character",   metavar="character",   default='100',  help="The number of gene modules to identify"),
  make_option(c("-t", "--knn_threshold"),         type = "character",   metavar="character",   default='0.8',  help="The minimum correlation value allowed to display connections between modules (not cells). Values should range between 0 and 1. Values around 0.3 to 0.6 reflect high correlation in single cells."),
  make_option(c("-a", "--assay"),                 type = "character",   metavar="character",   default='RNA',  help="Assay to be used in the analysis."),
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
suppressMessages(suppressWarnings({
  library(Seurat) 
  library(rafalib)
  library(WGCNA)
  library(circlize)
  library(RANN)
}))
#---------



#############################
### LOAD Seurat.v3 OBJECT ###
#############################
cat("\n### LOADING Seurat.v3 OBJECT ###\n")
DATA <- readRDS(opt$Seurat_object_path)
DATA@active.assay <- opt$assay
#---------



#####################
### PERFORM WGCNA ###
#####################
cat("\n### Performing WGCNA ###\n")

if(file.exists(paste0(opt$output_path,"/hclust_object.rds"))){
  h <- readRDS(paste0(opt$output_path,"/hclust_object.rds"))
} else {
  if( opt$gene_list_use == "none"){
    try(gene_list_use <- unique( object@assays[[assay]]@var.features ))
  } else {
    try(gene_list_use <- unique(rownames(read.csv2(opt$gene_list_use, row.names = 1))) )
  }
  gene_list_use <- gene_list_use[ gene_list_use %in% rownames(DATA@assays[[ opt$assay ]]@data) ]
  print(gene_list_use)
  #
  rowmax <- apply(DATA@assays[[ opt$assay ]]@data[ gene_list_use,],1,function(x) sum( x>0 ) )
  
  
  TOM <- WGCNA::cor(Matrix::t(DATA@assays[[ opt$assay ]]@data [ names(rowmax)[ rowmax >= 3 ] ,]), nThreads = 0)
  dim(TOM)
  TOM <- WGCNA::cor(1 - TOM, nThreads = 0)
  # TOM[is.na(TOM)] <- 0
  d <- as.dist( 1 - TOM )
  h <- hclust(d, method = "ward.D2")
  saveRDS(h,paste0(opt$output_path,"/hclust_object.rds"))
}
modules <- cutree(h, k = as.numeric(opt$number_of_modules) )
write.csv2(modules, paste0(opt$output_path,"/gene_modules_",opt$number_of_modules,".csv"),row.names = T)
#---------



############################
### COMPUTE MODULE MEANS ###
############################
cat("\n### Calculating gene module average expression ###\n")

rowmax <- apply(DATA@assays[[ opt$assay ]]@data[names(modules),],1,function(x) sort(x[x>0],T)[1] )
module_means <- rowsum(as.matrix(DATA@assays[[ opt$assay ]]@data[names(modules),]/(rowmax+1) ), modules)
module_means <- t(module_means)/c(table(modules))
DATA@assays[["modules"]] <- Seurat::CreateAssayObject(data = t(as.matrix(module_means)), min.cells = 0,min.features = 0)
#---------




########################
### MODULE KNN GRAPH ###
########################
cat("\n### Creating a gene module KNN-graph ###\n")


if(file.exists(paste0(opt$output_path,"/gene_modules_KNN_filt.csv"))){
  NN <- read.csv2( paste0(opt$output_path,"/gene_modules_KNN_filt.csv") , row.names = 1)
} else {
  NN <- WGCNA::cor( t(DATA@assays$modules@data) , nThreads = 0)
  # NN[is.na(NN)] <- 0
  NN <- RANN::nn2(NN, k = 5, eps = 0)
  NN <- data.frame( rep(NN$nn.idx[,1], ncol(NN$nn.idx)-1), c(NN$nn.idx[,-1]), c(NN$nn.dists[,-1])/2 )
  colnames(NN) <- c("from","to","dist")
  write.csv2(NN, paste0(opt$output_path,"/gene_modules_KNN_all.csv"),row.names = T)
  NN <- NN[ as.numeric(NN$dist) > as.numeric(opt$knn_threshold) , ]
  NN$weight <- (1 - NN$dist) / 2
  NN$scaled_weight <- ((1-NN$weight) - min(1-NN$weight) ) / (max(1-NN$weight) - min(1-NN$weight) )
  print(dim(NN))
  write.csv2(NN, paste0(opt$output_path,"/gene_modules_KNN_filt.csv"),row.names = T)

}
# dend <- hclust( dist( a@assays$modules@data ) ,method = "complete" )
# dend2 <- hclust( as.dist( 1-cor( t(a@assays$modules@data) )) ,method = "complete" )
# modules_interest <- sort(unique(c(32,23,28,2,3,13,60,4,26,70,7,6,50,17,65,41,11,23,28,34,54,25)))
if(opt$modules_of_interest != "none"){
  opt$modules_of_interest <- sort(unique(as.numeric(unlist(strsplit(opt$modules_of_interest,",")))))
}
#---------




#############################################
### DEFINING MODULE RELEVANCE TO METADATA ###
#############################################
cat("\n### Defining the relevance of each module ###\n")
n <- as.numeric(opt$number_of_modules)

if(opt$columns_metadata[1] != "none"){
  meta_data <- sort(unique(unlist(strsplit(opt$columns_metadata,","))))
  print(meta_data)

  for(x in meta_data){
    message(x)
    dev <- sapply(1:n, x=x, function(i,x){mod <- glm( DATA@assays$modules@data[i,] ~ DATA@meta.data[,x] ); c(mod$deviance, mod$null.deviance) } )
    assign( paste0("deltadeviance","_",x) , (dev[2,]-dev[1,])/dev[2,])
  }
  
  #PLOT per module
  pdf(paste0(opt$output_path,"/Module_relevance_plot.pdf"),width = 10,height = 2*length(meta_data),useDingbats = F)
  rafalib::mypar(length(meta_data),1,mar=c(2,10,1,1),mgp = c(2,0.5,0) )
  for(x in meta_data){
    i <- get(paste0("deltadeviance_",x))
    barplot(i,las=2, pch= 16,col="grey80",
            xaxs="i",names.arg = 1:n,border = NA,ylim=c(0,0.6),cex.main=.7)
    if(opt$modules_of_interest[1] != "none"){
    text( (opt$modules_of_interest*1.2-0.5), i[ opt$modules_of_interest ], labels = "*",pos = 3 )}
    mtext(x,side = 2,las=1,cex=.8,adj = 1,line = 2)
    lines( c(0,n*1.2), c(0,0) )
    lines( c(0,n*1.2), 0.1*c(1,1), lty=2, lwd=.5 )
    points( (1:n)*1.2-0.5 ,i, las=1, cex=.8, pch= 21,bg="grey70")
  }
  dev.off()
}


#---------




########################
### PLOT CIRCOS PLOT ###
########################
cat("\n### Plotting circos plots ###\n")

pdf(paste0(opt$output_path,"/Circos_plot.pdf"),6,6,useDingbats = F)

rafalib::mypar(mar=c(0,0,0,0))
circos.par(start.degree = 90)

circos.initialize("foo",xlim = c(0,n+10 ))
pos = circlize(x = 0:( n+7 ), 0:( n+7 ), sector.index = "foo")

# ADD HIGHLIGHT ON MODULES OF INTEREST WITH ASTERISK
if(opt$modules_of_interest != "none"){
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text((1:n)-0.25, rep(0, n), ifelse(1:n %in% opt$modules_of_interest,"*",""), col = "black",
  facing = "outside", adj = c(0, 0.5),cex = .9)},
  bg.border = NA, track.height = 0.05,track.margin=c(0.00,0.00),cell.padding=c(0,0,0,0))
}

# ADD IMPORTANCE OF EACH MODULE TO EACH METADATA
for(x in meta_data){
  message(x)
  
  circos.track(ylim = c(0, .5), panel.fun = function(z, y) {
  circos.rect(1:n-.1, rep(0, n), 1:n-.8, get(paste0("deltadeviance_",x)) ,
              col = colorRampPalette(c("grey95","black"))(100)[ round(get(paste0("deltadeviance_",x))*98)+1 ], border = NA)
              }, bg.border = NA, track.height=0.05,track.margin=c(0,.03),cell.padding=c(0,0,0,0))
  circos.yaxis(side="left",at = c(0.5),labels = paste0("delta_",x ), labels.cex = .7,col = "white")
}

# ADD THE NUMBER OF GENES IN EACH MODULE
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.rect(1:n-.1, rep(0, n), 1:n-.8, (table(modules)/max(table(modules))) , col = scales::hue_pal()(n), border = NA)}, 
  bg.border = NA, track.height=0.05,track.margin=c(0,.01),cell.padding=c(0,0,0,0))
circos.yaxis(side="left",at = c(0,1),labels = c(0,max(table(modules))), labels.cex = .7)
circos.yaxis(side="left",at = c(0.5),labels = "no. of genes        ", labels.cex = .7)


# ADD MODULE NAMES
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(1:n-0.5, rep(0, n), 1:n, col = scales::hue_pal()(n),
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = .7)
}, bg.border = NA, track.height = 0.05,track.margin=c(0,0.02),cell.padding=c(0,0,0,0))
circos.yaxis(side="left",at = c(0.5),labels = "module ID", labels.cex = .6,col = "white")


# ADD THE KNN GRAPH OF MODULE SIMILARITY
for(i in 1:nrow(NN)){
  circos.link("foo", NN[i,1]-0.5, 'foo', NN[i,2]-.5,h.ratio = 1,
              lwd = NN[i,4]*2+.2,col = paste0( scales::hue_pal()(n)[min(NN[i,1:2])],"60")  )}
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.05,track.margin=c(0.00,0.04),cell.padding=c(0,0,0,0))
circos.yaxis(side="left",at = c(0.5),labels = "module kNN", labels.cex = .7,col = "white")

circos.clear()
dev.off()
#---------








###############################
### SAVING Seurat.v3 OBJECT ###
###############################
cat("\n### Saving Seurat object ###\n")
saveRDS(DATA, file = opt$Seurat_object_path )
#---------



#############################
### SYSTEM & SESSION INFO ###
#############################
#---------
print_session_info()
#---------
