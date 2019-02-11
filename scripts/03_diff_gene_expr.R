#!/usr/bin/env Rscript

### LOAD LIBRARIES
#---------
library(optparse)
#---------


### DEFINE PATH TO LOCAL FILES
#---------
cat("\nRunning scQC with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object FILE."),
  make_option(c("-c", "--clustering_use"),        type = "character",   metavar="character",   default='none',  help="The clustering column to be used. Should be chosen from one of the method in script 02."),
  make_option(c("-m", "--metadata_use"),          type = "character",   metavar="character",   default='none',  help="A metadata column name can be supplied to be compared for differential expression for each clusters, as a second step. If nothing or incorrect values are provided, this step is ignored. Multible parameters can be provided and should be comma separated: 'TERM1,TERM2'. In this case, the analysis will be done for each factor individually."),
  make_option(c("-e", "--exclude_cluster"),       type = "character",   metavar="character",   default='none',  help="Clusters to be exluded from the analysis. Usefull for removing outlier cells. If the cluster name is not found, this will be ignored. Multible parameters can be provided and should be comma separated: 'TERM1,TERM2'. In this case, both clusters will be excluded from the DGE analysis."),
  make_option(c("-f", "--aux_functions_path"),    type = "character",   metavar="character",   default='none',  help="File with supplementary functions"),
  make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output DIRECTORY.")
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
pkgs <- c("rafalib","dplyr","RColorBrewer","scales","igraph")
inst_packages(pkgs)
#---------




### LOAD Seurat OBJECT 
#---------
cat("\nLoading/ data and metadata ...\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------



### Finding differentially expressed genes (cluster biomarkers)
#---------
DATA@ident <- factor(NULL)
DATA <- SetIdent(DATA,ident.use = DATA@meta.data[,opt$clustering_use])

#If the cluster to be excluded is present in the data, it will be removed
if(sum(as.character(unlist(strsplit(opt$exclude_cluster,","))) %in% unique(DATA@meta.data[,opt$clustering_use])) > 0 ){
  DATA <- SubsetData(DATA, cells.use = DATA@cell.names[! (DATA@meta.data[,opt$clustering_use] %in% as.character(unlist(strsplit(opt$exclude_cluster,","))) )]) #Filter out cell with no assigned clusters
  DATA@meta.data <- DATA@meta.data[DATA@cell.names[! (DATA@meta.data[,opt$clustering_use] %in% as.character(unlist(strsplit(opt$exclude_cluster,","))) )],]
}

if(file.exists(paste0(opt$output_path,"/Cluster_marker_genes.csv"))){
  cat("\nThe following Marker gene list was found here:\n",paste0(opt$output_path,"/Cluster_marker_genes.csv"),"\nDifferential expression among clusters will be skiped\n")
  DATA_markers <- droplevels(read.csv(paste0(opt$output_path,"/Cluster_marker_genes.csv"),row.names = 1))
} else {
  cat("\nThe Marker gene list was not found in the output folder. Computing differential expression among clusters ...")

  cat("\nThe following clusters (",opt$clustering_use,") will be used for analysis ...\n")
  print(unique(DATA@meta.data[,opt$clustering_use]))
  
  #DATA <- BuildSNN(DATA,reduction.type = "tsne",plot.SNN = F,k.param = 3,prune.SNN = .1)
  DATA_markers <- FindAllMarkers(object = DATA, only.pos = T)
  write.csv(DATA_markers,file = paste0(opt$output_path,"/Cluster_marker_genes.csv"),row.names = T)
}
#---------




### Plot a tSNE and SNNgraph-connected tSNE
#---------
cat("\nPlotting tSNE and SNN ...\n")
png(filename = paste0(opt$output_path,"/FINAL_tSNE_clustering.png"),width = 700,height = 600,res = 150)
TSNEPlot(object = DATA,group.by="ident")
invisible(dev.off())

png(filename = paste0(opt$output_path,"/tSNE_SNN_cluster_plot.png"),width = 600,height =650,res = 100)
net <- graph.adjacency(adjmatrix = as.matrix(DATA@snn), mode = "undirected", weighted = TRUE, diag = FALSE)
plot.igraph(x = net, layout = as.matrix(x = DATA@dr$tsne@cell.embeddings), edge.width = E(graph = net)$weight/2, vertex.label = NA,
            edge.color = colorRampPalette(c("grey90","black"))(50)[round(E(graph = net)$weight*49+1)],
            vertex.size = 3,vertex.frame.color=hue_pal(l=50, c=80)(length(levels(DATA@ident)))[DATA@ident],
            vertex.color = hue_pal()(length(levels(DATA@ident)))[DATA@ident])
invisible(dev.off())
#---------




### Plot a heatmap with the top genes for each cluster
#---------
cat("\nPlotting heatmap of Cluster Marker genes ...\n")
DATA_markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
#pdf(paste0(opt$output_path,"/Cluster_markers_heatmap.pdf"),width = 10,height = 10, useDingbats = F)
png(filename = paste0(opt$output_path,"/Cluster_markers_heatmap.png"),width = 900,height = 900,res = 150)
DoHeatmap(object = DATA, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE, use.scaled = F)
invisible(dev.off())
#---------



### Plot violin plots with the top genes for each cluster
#---------
 png(filename = paste0(opt$output_path,"/violinPlot_genes_per_cluster_scaled.png"),width = 200*30,height = 200*2*length(as.character(unique(top10$gene)))/10,res = 150)
 VlnPlot(object = DATA, features.plot = as.character(unique(top10$gene)),point.size.use = .1,nCol=10)
 invisible(dev.off())
 
 png(filename = paste0(opt$output_path,"/tSNE_plot_genes_per_cluster_scaled.png"),width = 200*30,height = 200*3*length(as.character(unique(top10$gene)))/10,res = 150)
 FeaturePlot(object = DATA, features.plot = as.character(unique(top10$gene)), cols.use = c("grey", "blue"), reduction.use = "tsne",pt.size = 1,nCol=10)
 invisible(dev.off())
#---------



if(sum(as.character(unlist(strsplit(opt$metadata_use,","))) %in% colnames(DATA@meta.data)) > 0 ){
  for(k in as.character(unlist(strsplit(opt$metadata_use,",")))){
    
### Plot the percentage of metadata in each cluster and vice-versa
#---------
out <- paste0(opt$output_path,"/Analysis_",k)
if(!dir.exists(out)){dir.create(out,recursive = T)}
cat("\nPlotting cell distribution across\t",k," ...\n")
my_meta_levels <- as.character(unique(DATA@meta.data[,k]))

proportion <- as.data.frame(lapply(levels(DATA@ident),function(x){c(unname(table(DATA@meta.data[DATA@ident==x,k]))) [1:length(my_meta_levels)]} ))
proportion[is.na(proportion)] <- 0
rownames(proportion) <- my_meta_levels
colnames(proportion) <- levels(DATA@ident)

sa <- cbind(stack(as.data.frame(proportion/rowSums(proportion))), rep(rownames(proportion),ncol(proportion)) )
colnames(sa) <- c("prop","clusters",k)

proportion2 <- t(as.matrix(t(proportion)/colSums(proportion)))
cl_order <- order(proportion2[my_meta_levels[1],],proportion2[my_meta_levels[2],],decreasing=T)
sa2 <- cbind(stack(as.data.frame(t(t(proportion[,cl_order])/colSums(proportion[,cl_order])) )), rep(rownames(proportion),ncol(proportion)) )
colnames(sa2) <- c("prop","clusters",k)

assign(k,k)

png(filename = paste0(out,"/barplot_overall.png"),width = 1300,height = 500,res = 150)
print(plot_grid(ggplot(data=sa2,mapping = aes(x=clusters, y=prop, fill=get(k)) ) + geom_bar(stat="identity"),ncol = 2,
          ggplot(data=sa,mapping = aes(x=get(k), y=prop, fill=clusters) ) + geom_bar(stat="identity") ))
invisible(dev.off())




### Identifying relevant markers across embrionic age for each population
#---------
cat("Calculating Differential gene expression among ",k,"   for each cluster ...\n")
marker_list <- list()
cluster_data <- list()
for(i in unique(DATA@ident)){
  cat(paste0("... Processing cell cluster #",i," ..."),"\n")
  temp <- SubsetData(DATA, cells.use = DATA@cell.names[DATA@ident == i]) #Select cell from a cluster
  temp@ident <- factor(NULL)
  temp <- SetIdent(temp,ident.use = temp@meta.data[,k])
  
  #check if all conditions have cells for differential expression
  if( length(table(temp@ident)) == length(my_meta_levels)  ){
    
    #check if the differential expression was previously calculated and saved as a .csv file
    if(file.exists(paste0(out,"/DEGs_in_cluster",i,".csv"))){
      
      cat("The following DEG list was found here:\n",paste0(out,"/DEGs_in_cluster",i,".csv"),"\nDifferential expression among clusters will be skiped\n")
      temp_markers <- data.frame(read.csv2(paste0(out,"/DEGs_in_cluster",i,".csv"),row.names = 1))
      
      } else {
        cat("\nThe following DEG list was not found in the output folder. Computing differential expression ...")
      
        temp_markers <- FindAllMarkers(object = temp, only.pos = T)  
        temp_markers <- temp_markers[(temp_markers$p_val < 0.01)&(temp_markers$avg_logFC > 0.3),]
        temp_markers <- temp_markers[order(temp_markers$p_val),]
        marker_list[[i]] <- temp_markers
        cluster_data[[i]] <- as.matrix(temp@data)[temp_markers$gene,]
        write.csv2(temp_markers,paste0(out,"/DEGs_in_cluster",i,".csv"),row.names = T)
        #if(nrow(temp_markers) > 20){top_temp <- temp_markers[1:20,]}
        
        }
      
      temp_markers %>% group_by(cluster) %>% top_n(5, avg_logFC) -> top_temp
      
      #Print the top 20 differentially expressed genes
      png(filename = paste0(out,"/DEGs_in_cluster",i,".png"),width = 200*5,height = 200*1.5*length(as.character(unique(top_temp$gene)))/4,res = 150)
      print(VlnPlot(object = temp, features.plot = as.character(unique(top_temp$gene)), point.size.use = .1)) #it does not work if you don't have the print command in front of it!
      dev.off()
      
    } else {
      cat(paste0("... Cell cluster #",i," has no cells in at least one of the conditions ..."),"\n")
    }
}}}
#---------    



cat("\n!!! Script executed Sucessfully !!!\n")


### System and session information
#---------
cat("\n\n\n\n... SYSTEM INFORMATION ...\n")
Sys.info()

cat("\n\n\n\n... SESSION INFORMATION ...\n")
sessionInfo()
#---------