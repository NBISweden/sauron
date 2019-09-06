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
cat("\nRunning DIFFERENTIAL EXPRESSION with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object FILE."),
  make_option(c("-i", "--VDJ_annotation_path"),   type = "character",   metavar="character",   default='none',  help="Path to the VDJ analysis output. Each dataset should be inside a folder. Inside each folder should be the file 'filtered_contig_annotations.csv' from cellranger, or similar."),
  make_option(c("-e", "--top_TCRs"),              type = "character",   metavar="character",   default='none',  help="How many TCRs to inlcude for differential expression analysis"),
  make_option(c("-p", "--paired_only"),           type = "character",   metavar="character",   default='FALSE', help="Logical. Whether to output only TCR that have a pair."),
  make_option(c("-o", "--only_coding_cdr3"),      type = "character",   metavar="character",   default='TRUE',   help="Logical. Whether to output only coding CDR3 sequences"),
  make_option(c("-o", "--same_scale"),            type = "character",   metavar="character",   default='TRUE',   help="Logical. Whether to use the same scale when plotting expression withing clusters. The scale is define as all cells."),
  make_option(c("-a", "--assay"),                 type = "character",   metavar="character",   default='RNA',   help="Assay to be used in the analysis."),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output DIRECTORY.")
) 
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)

pal <- scales::hue_pal()
mypal <- c(hue_pal()(9),
         RColorBrewer::brewer.pal(9,"Set1"),
         RColorBrewer::brewer.pal(8,"Set2"),
         RColorBrewer::brewer.pal(12,"Set3"),
         RColorBrewer::brewer.pal(9,"Pastel1"),
         RColorBrewer::brewer.pal(8,"Pastel2"),
         RColorBrewer::brewer.pal(8,"Accent") ) 
#---------



##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")
initial.options <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(sub("--file=","",initial.options[grep("--file=",initial.options)]))
source( paste0(script_path,"/inst_packages.R") )
pkgs <- c("rafalib","dplyr","RColorBrewer","scales","igraph","pheatmap","Seurat","venn")
inst_packages(pkgs)
#---------



##########################
### LOAD Seurat OBJECT ###
##########################
cat("\nLoading/ data and metadata ...\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------



###############################
### LOAD TCR CLONOTYPE DATA ###
###############################
cat("\nLoading data and metadata ...\n")
datasets <- list.dirs(opt$VDJ_annotation_path,recursive = F,full.names = F)

cat("\nThe following VDJ datasets were found: ...\n")
print(datasets)

VDJ <- data.frame()
for(i in datasets){
  temp <- read.csv(paste0(opt$VDJ_annotation_path,"/",i,"/filtered_contig_annotations.csv") )
  temp$barcode <- paste0(sub("-.*","",temp$barcode),"_",i)
  temp$dataset <- i
  VDJ <- rbind(VDJ,temp)
  #assign(i, temp)
}
rm(temp)
gc()
#---------



###############################
### PRE-PROCESSING VDJ DATA ###
###############################
#Merge VDJ names
VDJ$merged <- sapply(1:nrow(VDJ),function(x){ a <- as.character(unlist(VDJ[x,c("v_gene", "d_gene",  "j_gene" ,"c_gene")])); paste(a,collapse = "_")})

cat("\n The dimentions of your UNFILTERED VDJ table filtering table is:\n")
dim(VDJ)
res <- data.frame()
for(i in unique(VDJ$dataset) ){ res <- rbind(res,table( VDJ$chain [ VDJ$dataset == i] )) } 
res <- rbind( res, TOTAL=table( VDJ$chain ) )
rownames(res) <- c(datasets,"TOTAL")
colnames(res) <- names(table( VDJ$chain ))
cat("\n The UNFILTERED VDJ table contains this amount of chains ...\n")
print(res)
write.csv2(res,paste0(opt$output_path,"/VDJ_chains_per_dataset_UNFILTERED.csv"))



###########################################
# Filter VDJs that do not code anything ###
###########################################
cat("\n Filtering productive CDR3 sequences with high confidence ...\n")
if( casefold( opt$only_coding_cdr3 ) %in% c("true","yes") ){
  VDJ <- VDJ[VDJ$cdr3 != "None",]
  VDJ <- VDJ[VDJ$productive == "True",]
  VDJ <- VDJ[VDJ$raw_consensus_id != "None",]
  VDJ <- VDJ[VDJ$high_confidence != "False",]}


cat("\n The dimentions of your FILTERED VDJ table filtering table is:\n")
dim(VDJ)
res <- data.frame()
for(i in unique(VDJ$dataset) ){ res <- rbind(res,table( VDJ$chain [ VDJ$dataset == i] )) } 
res <- rbind( res, TOTAL=table( VDJ$chain ) )
rownames(res) <- c(datasets,"TOTAL")
colnames(res) <- names(table( VDJ$chain ))
cat("\n The FILTERED VDJ table contains this amount of chains ...\n")
print(res)
write.csv2(res,paste0(opt$output_path,"/VDJ_chains_per_dataset_FILTERED.csv"))



####################################
### Split each chain in a object ###
####################################
for(i in as.character(unique(VDJ$chain)) ){  assign(i , VDJ[ VDJ$chain==i ,])  }

pdf(paste0(opt$output_path,"/VDJ_Ig_detection_overlap.pdf"),4,4,useDingbats = F)
try(ll <- list(IGH=IGH$barcode, IGK=IGK$barcode, IGL=IGL$barcode ))
try(a<-attr(venn(ll, zcolor = pal(7)[1:3],opacity = .2), "intersections"))
dev.off()

pdf(paste0(opt$output_path,"/VDJ_TCR_detection_overlap.pdf"),4,4,useDingbats = F)
try(ll <- list(TRA=TRA$barcode,  TRB=TRB$barcode, TRD=TRD$barcode, TRG=TRG$barcode))
try(a<-attr(venn(ll,zcolor = pal(7)[4:7] ,opacity = .2), "intersections"))
dev.off()

pdf(paste0(opt$output_path,"/VDJ_TCRab_and_Ig_detection_overlap.pdf"),4,4,useDingbats = F)
try(ll <- list(TRA=TRA$barcode,  TRB=TRB$barcode, IGH=IGH$barcode, IGK=IGK$barcode, IGL=IGL$barcode) )
try(a<-attr(venn(ll,zcolor = pal(7)[1:5],opacity = .2), "intersections"))
dev.off()

VDJ_reads <- sapply( unique(VDJ$barcode) , function(x) { sum(VDJ[ VDJ$barcode == x, c("reads") ])  })
VDJ_umis <- sapply( unique(VDJ$barcode) , function(x) { sum(VDJ[ VDJ$barcode == x, c("umis") ])  })
VDJ_dataset <- sub(".*_","",unique(VDJ$barcode))

pdf(paste0(opt$output_path,"/VDJ_umis_vs_reads_per_cell.pdf"),5,5,useDingbats = F)
plot(log(VDJ_reads),log(VDJ_umis),cex=.1,main="VDJ_umis_vs_reads_per_cell",col = pal(length(datasets)) [as.numeric(as.factor(VDJ_dataset))] )
dev.off()



##############################################
### match each cdr3 per chain to each cell ###
##############################################
cat("\nDefining clonotypes based on the most abundant CDR3 sequences ...\n")
VDJ$IGH <- IGH$cdr3[ match(VDJ$barcode, IGH$barcode) ]
VDJ$IGL <- IGL$cdr3[ match(VDJ$barcode, IGL$barcode) ]
VDJ$IGK <- IGK$cdr3[ match(VDJ$barcode, IGK$barcode) ]
VDJ$TRA <- TRA$cdr3[ match(VDJ$barcode, TRA$barcode) ]
VDJ$TRB <- TRB$cdr3[ match(VDJ$barcode, TRB$barcode) ]
VDJ$TRD <- TRD$cdr3[ match(VDJ$barcode, TRD$barcode) ]
VDJ$TRG <- TRG$cdr3[ match(VDJ$barcode, TRG$barcode) ]
#head(VDJ[,c("barcode","chain","cdr3","umis","IGH","IGL","IGK","TRA","TRB")],20)
#tail(VDJ[,c("barcode","chain","cdr3","umis","IGH","IGL","IGK","TRA","TRB")],20)



#################################################
### Define clonotypes based on CDR3 sequences ###
#################################################
VDJ$TCRab_clonotype_seq <- paste0(VDJ$TRA,"_",VDJ$TRB)
VDJ$TCRgd_clonotype_seq <- paste0(VDJ$TRG,"_",VDJ$TRD)
VDJ$IGhl_clonotype_seq <- paste0(VDJ$IGH,"_",VDJ$IGL)
VDJ$IGhk_clonotype_seq <- paste0(VDJ$IGH,"_",VDJ$IGK)
write.csv2( VDJ , paste0(opt$output_path,"/VDJ_table.csv") )



######################################
### Filter clonotypes without pair ###
######################################
if( casefold(opt$paired_only) %in% c("true","yes") ){
  for(i in c("TCRab_clonotype_seq","TCRgd_clonotype_seq","IGhl_clonotype_seq","IGhk_clonotype_seq")){
    VDJ <-  VDJ[ !grepl("NA[_]|[_]NA",VDJ[,i]) | grepl("NA[_]NA",VDJ[,i]), ]
  }
}
cat("\n The FILTERED VDJ table contains this amount of chains ...\n")
write.csv2( VDJ , paste0(opt$output_path,"/VDJ_table_paired_only.csv") )
print(dim(VDJ))



##################################
### Update Chain table objects ###
##################################
for(i in as.character(unique(VDJ$chain))) {  assign(i , VDJ[ VDJ$chain==i ,])  }
#---------




#################################################
### Define levels of analysis to iterate over ###
#################################################
ls <- c("v_gene","d_gene","j_gene","cdr3","cdr3_nt","merged","TCRab_clonotype_seq","TCRgd_clonotype_seq","IGhl_clonotype_seq","IGhk_clonotype_seq")
ls <- sapply(ls, function(x) { length( levels( factor( VDJ[,x]) ))  })
ls <- names(ls)[ls >= 2]

for(j in ls ){
#for(j in c("cdr3","merged","clonseq","clon")){
  cat("\nPROCESSING ANALYSIS FOR: ",j," ...\n")
  output_path <- paste0(opt$output_path,"/",j)
  if(!dir.exists(output_path)){dir.create(output_path,recursive = T)}
  
  if( grepl("clonotype",j) ){ pars <- c("TRA|TRB","IGH|IGL","IGH|IGK")
  }else{ pars <- c("TRA","TRB","IGH","IGK","IGL") }
  
  for(k in pars){
#---------------



#####################################
### Computing relative abundances ###
#####################################
cat(" Computing relative abundances ...\n")
    
temp_VDJ <- VDJ[ grepl( k , VDJ$chain ), ]
#temp_VDJ <- temp_VDJ[!duplicated(temp_VDJ$barcode),]
temp_VDJ <- temp_VDJ[temp_VDJ$barcode %in% rownames(DATA@meta.data),]
abund_all_cell <- sapply( unique(temp_VDJ$barcode) , function(p) {  temp_VDJ[temp_VDJ$barcode==p,j] [1] } )
abund_all <- sort(table( as.character(na.omit( abund_all_cell ))),T)
  
tot_all <- sum(abund_all)

write.csv2(cbind(counts=abund_all, percentage=100*abund_all/tot_all)[abund_all>0,],paste0(output_path,"/",j,"_abundance_TOTAL.csv"))
x_top <- sort(abund_all[1:min(31,sum(abund_all>0))],T)
x_top <- c( x_top,others = sum(abund_all[ ! names(abund_all) %in% names(x_top)] ) )
#---------------


#####################################
### Plotting TOTAL CDR3 abundance ###
#####################################
png(filename = paste0(output_path,"/",j,"_diversity_plots.png"),width = 3000, height = 1200,res = 150)
par(mfrow=c(1,2))
par(mar=c(10,18,10,12),xpd=F)

mycol <- c(mypal[1:(length(x_top)-1)],"grey90")

barplot(rev(x_top),horiz=T,las=1,cex.names=.8,yaxs="i",xaxs="i",xlim=c(0,max(x_top)*1.2),ann=FALSE,axes=FALSE,
        main=paste0(j," diversity"),xlab="number of cells",border=NA,col=rev(mycol),xpd=F)
lines(c(0,0),c(0,length(x_top)*1.2),xpd=T,lwd=1)
axis(1,at=seq(0,ceiling(max(x_top)/12)*12,length.out = 5),las=2)
axis(3,at=seq(0,ceiling(max(x_top)/10)*10,length.out = 5),las=2,labels = paste0(round(seq(0,ceiling(max(x_top)/12)*12/tot_all*100,length.out = 5),1),"%") )

par(mar=c(2,5,2,5))
pie(x_top,clockwise = T,main = paste0("\n\n",j," diversity"),col = mycol,labels = names(x_top)[1:9],radius = .6,border = NA,xpd=T,line=0)

invisible(dev.off())
#---------------




############################################################################
### Plotting METADATA-specific CDR3 abundance of TCRa and TCRb sequences ###
############################################################################
cat("   Plotting Cluster-specific CDR3 abundance of TCRa and TCRb sequences ...\n")
for(clustering in opt$metadata_use){
    for(i in sort(unique(DATA@meta.data[,clustering])) ){
    
    cell_use <- rownames(DATA@meta.data)[DATA@meta.data[,clustering] == i]
    clus_VDJ <- temp_VDJ[temp_VDJ$barcode %in% cell_use ,]

    abund_clus <- sapply( unique(clus_VDJ$barcode) , function(p) {  clus_VDJ[clus_VDJ$barcode==p , j] [1] } )
    abund_clus <- table( as.character(na.omit( abund_clus )) )
    
    if( casefold(opt$same_scale) %in% c("true","yes") ){
      x <- abund_all
      x[x>=0] <- 0
      x[match(names(abund_clus), names(x))] <- abund_clus
      abund_clus <- x
    } else { tot_all <- sum(abund_clus) }
    write.csv2(cbind(counts=abund_clus, percentage=100*abund_clus/sum(DATA@meta.data[,clustering] == i))[abund_clus>0,],paste0(output_path,"/",k,"_abundance_",clustering,"_",i,".csv"))
    clus_top <- abund_clus[1:min(31,sum(abund_clus>0))]
    clus_top <- c( clus_top, others = sum(abund_clus[ ! names(abund_clus) %in% names(x_top)] ) )
    
    
    png(filename = paste0(output_path,"/",j,"_diversity_plots_",clustering,"_",i,".png"),width = 3000, height = 1200,res = 150)
    par(mfrow=c(1,2))
    par(mar=c(10,18,10,12),xpd=F)
    mycol <- c(mypal[1:(length(clus_top)-1)],"grey90")
    
    barplot(rev(clus_top),horiz=T,las=1,cex.names=.8,yaxs="i",xaxs="i",xlim=c(0,max(x_top)*1.2),ann=FALSE,axes=FALSE,
            main=paste0(j," diversity"),xlab="number of cells",border=NA,col=rev(mycol),xpd=F)
    lines(c(0,0),c(0,length(x_top)*1.2),xpd=T,lwd=1)
    axis(1,at=seq(0,ceiling(max(x_top)/12)*12,length.out = 5),las=2)
    axis(3,at=seq(0,ceiling(max(x_top)/10)*10,length.out = 5),las=2,labels = paste0(round(seq(0,ceiling(max(x_top)/12)*12/tot_all*100,length.out = 5),1),"%") )
    
    par(mar=c(2,5,2,5))
    pie(x_top,clockwise = T,main = paste0("\n\n",j," diversity"),col = mycol,labels = names(x_top)[1:9],radius = .6,border = NA,xpd=T,line=0)
    
    invisible(dev.off())
  }
}    
#---------------



######################################
### Clone abundance mapped to tSNE ###
######################################
cat("   Mapping clone abundance mapped to tSNE ...\n")
calc_abund <- function(x,i){
  rest_clon <- names(sort(table(x[,j]),T)) [ sort(table(x[,j]),T) %in% i]
  return(as.character(x[x[,j] %in% rest_clon,"barcode"]))}
cell_names <- colnames(DATA)


#Clone abundance TCR-a
rest_clon <- cell_names %in% calc_abund(temp_VDJ,1:2)
rest_clon <- ifelse(cell_names %in% calc_abund(temp_VDJ,1),"1",
             ifelse(cell_names %in% calc_abund(temp_VDJ,2:3),"2-3",
             ifelse(cell_names %in% calc_abund(temp_VDJ,3:20000),">3",NA) ) )


DATA@meta.data[[paste0(k,"_",j)]] <- abund_all_cell[ match( rownames(DATA@meta.data), names(abund_all_cell) ) ]
DATA@meta.data[[paste0(k,"_",j,"_abundance")]] <-  as.numeric(abund_all [ match( DATA@meta.data[[paste0(k,"_",j)]], names(abund_all) )])
  
  
DATA@meta.data[[paste0(k,"_",j,"_abundance")]] <- rest_clon
cat("   summary statistics ...\n")
print(table(rest_clon)) ; print(sum(is.na(rest_clon)))

if( !is.null( names(DATA@reductions)) ){
  for(red in names(DATA@reductions)){
    png(filename = paste0(output_path,"/",k,"_",j,"_abundance_",red,".png"),width = 900,height = 800,res = 150)
    print( FeaturePlot(DATA,reduction = red,cols = c("grey90","navy"),features = paste0(k,"_",j,"_abundance")) )
    dev.off()
    
    png(filename = paste0(output_path,"/",k,"_",j,"_abundance2_",red,".png"),width = 1800,height = 800,res = 150)
    print( DimPlot(DATA,reduction = red, group.by = paste0(k,"_",j), plot.title=paste0(k,"_",j)) )
    dev.off()
}}
#-----------------



#############################################
### Mapping TCR sequences onto tSNE plot ###
#############################################
if( !is.null( names(DATA@reductions)) ){
  for(red in names(DATA@reductions)){
    png(filename = paste0(output_path,"/",k,"_",j,"_ind_abund_",red,".png"),width = 600*5, height = 400*ceiling(length(x_top)/5),res = 150)
    par(mar=c(1,1,1,10),mfrow=c(ceiling(length(x_top)/5),5))
    for(i in names(x_top) ){
      clon_A <- ifelse( DATA@meta.data[[paste0(k,"_",j)]] == i , i , "Other")
      clon_A[!(cell_names %in% get(k)$barcode)] <- "NA"

      o <- order(factor(clon_A,levels = c(i,"Other","NA")),decreasing = T)
      
      plot(DATA@reductions[[red]]@cell.embeddings[o,],col=c("red","grey80","grey50")[factor(clon_A[o],levels = c(i,"Other","NA"))],pch=16,cex=.7,line=0,axes=F,main=paste0(i))
      legend(max(DATA@reductions[[red]]@cell.embeddings[o,1]),max(DATA@reductions[[red]]@cell.embeddings[o,2])*1.2,y.intersp = .9,cex = .95,
             legend = c(i,"Other","NA"),bty = "n",col =c("red","grey80","grey50"),pch = 16,pt.cex = 1.5,xpd=T)
    }
    invisible(dev.off())
}}
#---------------



# MSLLTEVETPTRNEWECRCSDSS
# 
# sort(table(TCR[TCR$cdr3=="CASSENSGNTLYF","cdr3_nt"]),T)[1:4]
# a <- cell_names %in% TCR[TCR$cdr3=="CASSENSGNTLYF","cdr3_nt"]
# table(DATA@meta.data$days_post_infection[a])
# 
# 
# a <- cell_names %in% TCR[TCR$cdr3_nt=="TGTGCCAGCAGTGAAAATTCTGGAAATACGCTCTATTTT","barcode"]
# table(DATA@meta.data$days_post_infection[a])
# 
# 
# a <- cell_names %in% TCR[TCR$cdr3_nt=="TGTGCCAGCAGTGAGAATTCTGGAAATACGCTCTATTTT","barcode"]
# table(DATA@meta.data$days_post_infection[a])
# 



#########################################
### TCR sequence co-detection heatmap ###
#########################################
if( !( j %in% c("clonseq","clon")) ){
cat("   Computing TCR sequence co-detection ...\n")
for(k in c("TCRA","TCRB")){
  temp <- names(sort(table(get(k)[,j]),T)[1:31])
  tcrs <- names(sort(table(get(k)$cdr3),T)[1:31])
  res <- data.frame(row.names = temp)
  for(i in tcrs){
    sel <- as.character(get(k)$barcode[get(k)$cdr3 == i])
    sel <- TCR[(TCR$barcode %in% sel) & (TCR$chain == ifelse(k=="TCRB","TRB","TRA")),]
    sel <- sel[sel[,j] %in% temp,]
    sel <- table(factor(sel[,j],levels = temp) )[temp]
    res <- cbind(res,as.numeric(sel) )
  }
  colnames(res) <- tcrs
  
  #plotting heatmap
  x <- pheatmap(res,fontsize_row = 8,fontsize_col = 8,border="grey90",main = paste0(k,"-to-",j," co-detection"),
           col=c("grey95",colorRampPalette(c("grey80","orange3","firebrick","red"))(400)),cellwidth = 8,cellheight = 8,angle_col = 90,
           filename = paste0(output_path,"/",k,"_cdr3-to-",j,"_co-detection_heatmap.png") )
  write.csv2(res[x$tree_row$order,x$tree_col$order],paste0(output_path,"/",k,"_cdr3-to-",j,"_co-detection.csv"))
  }
}
#---------------



##############################
### TCR chain co-detection ###
##############################
if( !( j %in% c("clonseq","clon")) ){
cat("   Computing TCR chain co-detection ...\n")

temp <- names(sort(table(TCRA[,j]),T)[1:31])
tcrs <- names(sort(table(TCRB[,j]),T)[1:31])
res <- data.frame(row.names = temp)
for(i in tcrs){
  sel <- as.character(TCRB$barcode[TCRB[,j] == i])
  sel <- TCR[(TCR$barcode %in% sel) & (TCR$chain == "TRA"),]
  sel <- sel[sel[,j] %in% temp,]
  sel <- table(factor(sel[,j],levels = temp) )[temp]
  res <- cbind(res,as.numeric(sel) )
}
colnames(res) <- tcrs
res <- res / sum(res)*100

#plotting heatmap
x <- pheatmap(res,fontsize_row = 8,fontsize_col = 8,border="grey90",main = paste0(k,"-to-",j," co-detection"),clustering_method = "ward.D2",
         col=c("grey95",colorRampPalette(c("grey80","orange3","firebrick","red"))(400)),cellwidth = 8,cellheight = 8,angle_col = 90,
         legend_breaks = seq(0,100,length.out = 11), filename = paste0(output_path,"/TCRA-TCRB_",j,"_co-detection_heatmap.png"))
write.csv2(res[x$tree_row$order,x$tree_col$order],paste0(output_path,"/TCR_chains_",j,"_co-detection.csv"))
}
#---------------



############################################
### TCR sequence to cluster co-detection ###
############################################
cat("   Computing TCR sequence to cluster co-detection ...\n")
for(clustering in opt$metadata_use){
  for(k in c("TCRp","TCRA","TCRB")[ifelse(j %in% c("clonseq","clon"),1,-1)]){
    tcrs <- names(sort(table(get(k)[,j]),T))
    res <- data.frame(row.names = tcrs)
    for(i in sort(unique(DATA@meta.data[,clustering]))){
      sel <- as.character(cell_names[DATA@meta.data[,clustering] == i])
      sel <- TCR[(TCR$barcode %in% sel) & (TCR$chain %in% ifelse(k == "TCRB","TRB",ifelse(k == "TCRA","TRA",c("TRB","TRA")))),]
      sel <- sel[sel[,j] %in% tcrs,]
      sel <- table(factor(sel[,j],levels = tcrs) )[tcrs]
      res <- cbind(res,as.numeric(sel) )
    }
    colnames(res) <- sort(unique(DATA@meta.data[,clustering]))
    res <- t(t(res) / colSums(res))*100
    
    pheatmap(res[1:31,],fontsize_row = 8,fontsize_col = 8,border="grey90",main = paste0("Cluster-",k," co-detection (% per cluster)"),cluster_cols = F,
             col=c("grey95",colorRampPalette(c("grey80","orange3","firebrick","red"))(400)),legend_breaks = seq(0,100,length.out = 11), cellwidth = 8,cellheight = 8,angle_col = 90,filename = paste0(output_path,"/Clusters-to-",k,"_co-detection_heatmap.png") )
  }
}
#---------------



##############################################################################
### Find differentially expressed genes for CASSENSGNTLYF vs CASSDNSGNTLYF ###
##############################################################################
cat("   Computing differentially expressed genes for CASSENSGNTLYF vs CASSDNSGNTLYF ...\n")
n_tcr <- 1:as.numeric(opt$top_TCRs)
n_tcr <- c("CASSENSGNTLYF","CASSDNSGNTLYF")

for(k in c("TCRp","TCRA","TCRB")[ifelse(j %in% c("clonseq","clon"),1,-1)]){
  tcrs <- names(sort(table(get(k)[,j]),T)[n_tcr])
  
  tcr_annot <- rep("Other",length(cell_names))
  for(i in tcrs){
    tcr_annot[cell_names %in% as.character(get(k)[  get(k)[,j] == i , "barcode"]) ] <- i }
  
  if(casefold(opt$include_single) == "true"){
  tcr_annot[cell_names %in% as.character(get(k)[  get(k)[,j] %in% names(table(get(k)[,j])[table(get(k)[,j]) == 1]) , "barcode"]) ] <- "single"}
  DATA@meta.data[,paste0(k,"_",j)] <- tcr_annot
  DATA@ident <- factor(NULL)
  DATA <- SetIdent(DATA,ident.use = tcr_annot )
  DATA@ident <- factor(DATA@ident,levels = c(tcrs,"Other","single") )
  temp <- SubsetData(DATA, cells.use = DATA@cell.names[DATA@ident != "Other"]) #Select cell from a cluster
  print(table(temp@ident))
  
  diff_exp <- FindAllMarkers(temp,print.bar = T,min.diff.pct = .1,min.pct = .2)
  write.csv2(diff_exp,file = paste0(output_path,"/",k,"_marker_genes.csv"),row.names = T)
  diff_exp %>% group_by(cluster) %>% top_n(20, avg_logFC) -> top10
  
  png(filename = paste0(output_path,"/violinPlot_",k,"_enriched_genes_CASSENSGNTLYF_CASSDNSGNTLYF.png"),width = 500*(2.5*length(sort(unique(temp@meta.data[,paste0(k,"_",j)])))+1),height = 500*ceiling((1+length(unique(top10$gene)))/10),res = 300)
  par(mfrow=c(ceiling((1+length(unique(top10$gene)))/10),10))
  for(l in unique(top10$gene)){
    violins(data = temp,gene = l,clustering = paste0(k,"_",j),plot_points = T,plot_axis = F)
  }
  par(mar=c(0,0,0,0))
  plot(0,0,type="n",yaxt="n",xaxt="n",frame.plot = F,ylab="")
  legend(0,0,y.intersp = .9,cex = .90,xjust = .5,yjust = .5,
         legend = sort(unique(temp@meta.data[,paste0(k,"_",j)])),bty = "n",col =hue_pal()(length(sort(unique(temp@meta.data[,paste0(k,"_",j)])))),pch = 16,pt.cex = 1.5,xpd=T)
  #VlnPlot(object = temp, features.plot = as.character(unique(top10$gene)),point.size.use = .1,nCol=10)
  invisible(dev.off())
}
#---------------



########################################################################################
### Find differentially expressed genes same TCR sequence between different clusters ###
########################################################################################
cat("   Computing differentially expressed genes same TCR sequence between different clusters ...\n")
for(clustering in opt$metadata_use){
  if(!dir.exists(paste0(output_path,"/",clustering))){dir.create(paste0(output_path,"/",clustering))}
  n_tcr <- 1:as.numeric(opt$top_TCRs)
  
  for(k in c("TCRp","TCRA","TCRB")[ifelse(j %in% c("clonseq","clon"),1,-1)]){
    tcrs <- names(sort(table(get(k)[,j]),T)[n_tcr])
    
    tcr_annot <- rep("Other",length(cell_names))
    for(i in tcrs){
      tcr_annot[cell_names %in% as.character(get(k)[  get(k)[,j] == i , "barcode"]) ] <- i }
    
    if(casefold(opt$include_single) == "true"){
      tcr_annot[cell_names %in% as.character(get(k)[  get(k)[,j] %in% names(table(get(k)[,j])[table(get(k)[,j]) == 1]) , "barcode"]) ] <- "single"}
    
    DATA@meta.data[,paste0(k,"_",j)] <- tcr_annot
    
    for(m in tcrs){
      cat(m)
      temp <- SubsetData(DATA, cells.use = DATA@cell.names[tcr_annot == m]) #Select cell from a cluster
      temp@ident <- factor(NULL)
      temp <- SetIdent(temp,ident.use = temp@meta.data[,clustering] )
      print(table(temp@ident))
      
      if( file.exists(paste0(output_path,"/",clustering,"/",k,"_marker_genes_",m,"_",clustering,".csv")) ){
        diff_exp <- read.csv2(paste0(output_path,"/",clustering,"/",k,"_marker_genes_",m,"_",clustering,".csv"),row.names = 1)
      } else {
        diff_exp <- FindAllMarkers(temp,print.bar = T,min.diff.pct = .1,min.pct = .2)}
      
      if(dim(diff_exp)[1]!=0){
        write.csv2(diff_exp,file = paste0(output_path,"/",clustering,"/",k,"_marker_genes_",m,"_",clustering,".csv"),row.names = T) 
        diff_exp %>% group_by(cluster) %>% top_n(20, avg_logFC) -> top10
    
        png(filename = paste0(output_path,"/",clustering,"/violinPlot_",m,"_DEGs_between_",clustering,".png"),width = 500*(2.5*length(sort(unique(temp@meta.data[,clustering])))+2),height = 500*ceiling((1+length(unique(top10$gene)))/10),res = 300)
        par(mfrow=c(ceiling((1+length(unique(top10$gene)))/10),10))
        for(l in unique(top10$gene)){
          violins(data = temp, gene = l, clustering = clustering, plot_points = T, plot_axis = T,col =hue_pal()(length(table(DATA@ident)))[ as.numeric(names(table(temp@ident)))+1] )
        }
        invisible(dev.off())
      } else {cat("  No differentially expressed genes found ...\n")}
    }
  }
}
#---------------
}



#Find differentially expressed genes for a particular TCR
#---------------
i <- "CASSENSGNTLYF"
DATA@ident <- factor(NULL)
id <- DATA@meta.data$TCRb_CASSENSGNTLYF
id[id == "NA"] <- "Other"
data <- SetIdent(data,ident.use = id )
DATA@meta.data$ident <- id
diff_exp <- FindMarkers(data,ident.1 = i, ident.2 = NULL,print.bar = T,min.diff.pct = .1,min.pct = .2)


png(filename = paste0(opt$output_path,"/DiffExpression_for_",i,".png"),width = 200*15,height = 200*2*ceiling(length(rownames(diff_exp)[1:20])/5),res = 200)
#VlnPlot(data,features.plot = rownames(diff_exp)[1:20],nCol = 5,point.size.use = .001)
par(mfrow=c(ceiling(length(rownames(diff_exp)[1:20])/5),5))
for(j in rownames(diff_exp)[1:20]){
  violins(data,j,"ident",T,F)
}
invisible(dev.off())
#---------------



#Plot the most abundant pair for certain TCR
#---------------
i <-"CASSENSGNTLYF"
  
sel <- TCR$barcode[TCR$cdr3 == i]
sel <- TCR[(TCR$barcode %in% sel) & (TCR$chain != "TRB"),]
sel <- sel[sel$cdr3 != "None",]
sel <- sort(table(as.character(sel$cdr3) ),T)

par(mar=c(4,8,2,2),mfrow=c(1,1))
barplot(sort(sel,F),horiz=T,las=1,cex.names=.65,yaxs="i",xaxs="i",cex.main=.8,
        border=NA,main=paste("TCRa CDR3 sequences\n to TCRb_",i),xlab="frequency")
#---------------






###################################
### SAVING RAW Seurat.v3 OBJECT ###
###################################
cat("\nSaving the RAW Seurat object ...\n")
saveRDS(DATA, file = paste0(opt$output_path,"/Seurat_object_TCR.rds") )
#---------



#############################
### SYSTEM & SESSION INFO ###
#############################
cat("\n\n\n\n... SYSTEM INFORMATION ...\n")
Sys.info()

cat("\n\n\n\n... SESSION INFORMATION ...\n")
sessionInfo()
#---------




# 
# ifelse()
# abundances <- cell_names
# abundances <- cell_names


# 
# 
# data@meta.data$tcrb_cdr3 <- TCRB[ match(cell_names,TCRB$barcode),"cdr3"]
# data@meta.data$tcrb_cdr3[is.na(data@meta.data$tcrb_cdr3)] <- "None"
# 
# sort(table(data@meta.data$tcrb_cdr3),T)[1:10]
# data@meta.data$clon_A <- ifelse(as.character(data@meta.data$tcrb_cdr3) == "CASSENSGNTLYF","clon_A","None")
# TSNEPlot(data,group.by="clon_A")
# 
# 
# cdr3 <- setNames( TCR[ match(cell_names,TCR$barcode),"raw_clonotype_id"],rownames(data@meta.data))
# clonotype <- setNames( TCR[ match(cell_names,TCR$barcode),"raw_clonotype_id"],rownames(data@meta.data))
# data <- AddMetaData(object = data, metadata = clonotype, col.name = "clonotype")
# data@meta.data[,"clonotype"]
# 
# 
# data@meta.data$clon1 <- ifelse(as.character(data@meta.data$clonotype) == "clonotype1","clon1","NA")
#   
# table(data@meta.data$clon1)
# sort(table(TCR$raw_clonotype_id),T)[1:10]
# 
#   
# TSNEPlot(data,group.by="clon1")
# 
# 


# 
# sort(table(clonotype),T)
# 



