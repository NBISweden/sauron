#!/usr/bin/env Rscript
### ANALYSIS OF FIBROBLASTS FOR M_KASPER_234 PROJECT
### PAULO CZARNEWSKI


### LOAD LIBRARIES
#---------
source("/Users/paulo.barenco/Desktop/NBIS/Scripts/czarlib/inst_packages.R")
pkglist <- c("Seurat","ggplot2","dplyr","biomaRt","scran","devtools","sva","RColorBrewer","tsne","Rtsne","flowPeaks","dbscan","vegan")
inst_packages(pkglist)

library(Seurat)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(scran)
library(sva)
library(RColorBrewer)
library(devtools)
library(scales)
#---------


### DEFINE PATH TO LOCAL FILES
#---------
Seurat_object_path <- "/Users/paulo.barenco/Desktop/NBIS/Projects/LIN/3-Clustering/Clustered_Seurat_object.rds"
clustering_level <- "dbscan.8"
clusters_to_select <- as.character( c(1) )

output_path <- "/Users/paulo.barenco/Desktop/NBIS/Projects/LIN/4-Analysis_Tcells"
if(!dir.exists(output_path)){dir.create(output_path)}
setwd(output_path)

col_scale <- c("grey85","navy")
if(!dir.exists(paste0(output_path,"/tSNE_plots"))){dir.create(paste0(output_path,"/tSNE_plots"))}
if(!dir.exists(paste0(output_path,"/PCA_plots"))){dir.create(paste0(output_path,"/PCA_plots"))}

#---------




### LOAD Seurat OBJECT (already filtered for genes with lowly-expressed genes and low-quality cells)
#---------
DATA <- readRDS(Seurat_object_path)
#---------




### Selecting the cells/clusters of interest
#---------
cells_use <- rownames(DATA@meta.data)[DATA@meta.data[,clustering_level] %in% clusters_to_select]

DATA2 <- CreateSeuratObject(DATA@raw.data[,cells_use], meta.data = DATA@meta.data[cells_use,])
DATA2 <- NormalizeData(DATA2)
DATA2 <- FindVariableGenes(DATA2,x.low.cutoff = 0.01,x.high.cutoff = 2,y.cutoff = 1,do.plot = F)
plot(log2(DATA2@hvg.info[,1]),DATA2@hvg.info[,3],cex=.1,col=ifelse(rownames(DATA2@hvg.info)%in%DATA2@var.genes,"red","black"),
     xlab="mean",ylab="Scaled dispersion")

vars.to.regress <- c("nUMI","SampleID","percent.mito","batch","S.Score","G2M.Score","sequencing_run","percent.Rpl","percent.Rps")
DATA2 <- ScaleData(object = DATA2, vars.to.regress = vars.to.regress)
#---------


### Running dimentionality reduction on the selected cells
#---------
DATA2 <- RunPCA(object = DATA2,pc.genes = DATA2@var.genes[1:300])

pdf(paste0(output_path,"/PCA_plots/PCA_leading_genes.pdf"),width = 10,height = 20,useDingbats = F)
VizPCA(object = DATA2, pcs.use = 1:20, nCol = 4,font.size = 1)
dev.off()

pdf(paste0(output_path,"/PCA_plots/PCA_plots.pdf"),width = 10,height = 9,useDingbats = F)
for(j in c("SampleID","batch","sequencing_run","days_post_infection","clean_up")){
  PCAPlot(object = DATA2, dim.1 = 1, dim.2 = 2, group.by=j,pt.size = .3) }
for(j in c("nUMI","nGene","S.Score","G2M.Score","percent.Rpl","percent.Rps")){
  FeaturePlot(object = DATA2, features.plot = j,cols.use = col_scale,reduction.use = "pca",pt.size = .3) }
dev.off()

DATA2 <- ProjectPCA(object = DATA2, do.print = T)

pdf(paste0(output_path,"/PCA_plots/PCA_genes_heatmap.pdf"),width = 10,height = 20,useDingbats = F)
PCHeatmap(object = DATA2, pc.use = 1:20, cells.use = 1000, do.balanced = TRUE, label.columns = FALSE)
dev.off()
#---------



### Run Non-linear dimensional reduction (tSNE)
#---------
DATA2 <- RunTSNE(object = DATA2, perplexity=30, max_iter=2000, theta=0.05, eta=1000,exaggeration_factor=12,dims.use = 1:8)

for(j in c("SampleID","batch","sequencing_run","days_post_infection","clean_up")){
  png(filename = paste0(output_path,"/tSNE_plots/tSNE_",j,".png"),width = 700,height = 600,res = 150)
  TSNEPlot(object = DATA2,group.by=j,pt.size = 1)
  dev.off()
}

for(j in c("nUMI","nGene","S.Score","G2M.Score","percent.Rpl","percent.Rps")){
  png(filename = paste0(output_path,"/tSNE_plots/tSNE_",j,".png"),width = 600,height = 600,res = 150)
  FeaturePlot(object = DATA2, features.plot = j, cols.use = col_scale,pt.size = 1)
  dev.off()
}

png(filename = paste0(output_path,"/tSNE_plots/tSNE_density.png"),width = 600,height =650,res = 100)
plot(DATA2@dr$tsne@cell.embeddings,pch=16,cex=0.7, col=densCols(DATA2@dr$tsne@cell.embeddings, colramp=colorRampPalette(c("grey80","grey80","navy", "dark green","orange", "red", "firebrick")), nbin = 500,bandwidth =7),las=1,xlab="tSNE_1",ylab="tSNE_2")
dev.off()

png(filename = paste0(output_path,"/tSNE_plots/tSNE_SNN_plot.png"),width = 600,height =650,res = 100)
DATA2 <- BuildSNN(DATA2,reduction.type = "tsne",plot.SNN = T,k.param = 3,prune.SNN = .1)
dev.off()
#---------




### Finding clusters using SNN
#---------
if(!dir.exists(paste0(output_path,"/clustering_SNN"))){dir.create(paste0(output_path,"/clustering_SNN"))}
for(k in seq(.1,2,by=.1)){
  if(k!=.1){
  DATA2 <- FindClusters(object = DATA2, reduction.type = "pca", dims.use = 1:8, resolution = k, print.output = F)
  } else { DATA2@ident <- factor(NULL) ; DATA2 <- FindClusters(object = DATA2, reduction.type = "pca", dims.use = 1:8, resolution = k, print.output = F, save.SNN = TRUE,force.recalc = T)}
  png(filename = paste0(output_path,"/clustering_SNN/tSNE_SNN_res=",format(k,nsmall=2),".png"),width = 700,height = 600,res = 150)
  TSNEPlot(object = DATA2, group.by=paste0("res.",k), pt.size = 1, plot.title= paste0("SNN clustering with res=",format(k,nsmall=2)))
  dev.off()
}
#---------



### Clustering using HDBSCAN
#---------
if(!dir.exists(paste0(output_path,"/clustering_HDBSCAN"))){dir.create(paste0(output_path,"/clustering_HDBSCAN"))}
for(k in seq(4,60,by=2)){
  clusters <- hdbscan(DATA2@dr$tsne@cell.embeddings,minPts = k)
  names(clusters$cluster) <- rownames(DATA2@meta.data)
  DATA2 <- AddMetaData(object = DATA2, metadata = clusters$cluster, col.name = paste0("hdbscan.",k))
  png(filename = paste0(output_path,"/clustering_HDBSCAN/tSNE_HDBSCAN_k=",format(k,nsmall=2),".png"),width = 700,height = 600,res = 150)
  TSNEPlot(object = DATA2, group.by=paste0("hdbscan.",k), pt.size = 1, plot.title= paste0("HDBSCAN clustering with minPts=",format(k,nsmall=2)))
  dev.off()
}
#---------



### Clustering using FLOWPEAKS
#---------
if(!dir.exists(paste0(output_path,"/clustering_FlowPeaks"))){dir.create(paste0(output_path,"/clustering_FlowPeaks"))}
for(k in seq(.1,2,by=.05)){
  c <- flowPeaks(DATA2@dr$tsne@cell.embeddings,h0 = k)
  cluster_ID <- assign.flowPeaks(c,DATA2@dr$tsne@cell.embeddings,tol = 0.02,fc = .1)
  names(cluster_ID) <- rownames(DATA2@meta.data)
  DATA2 <- AddMetaData(object = DATA2, metadata = cluster_ID, col.name = paste0("flowpeaks.",k))
  png(filename = paste0(output_path,"/clustering_FlowPeaks/tSNE_FlowPeaks_k=",format(k,nsmall=2),".png"),width = 700,height = 600,res = 150)
  TSNEPlot(object = DATA2, group.by=paste0("flowpeaks.",k), pt.size = 1, plot.title= paste0("FlowPeaks clustering with h0=",format(k,nsmall=2)))
  dev.off()
}
#---------



### Clustering using DBSCAN
#---------
if(!dir.exists(paste0(output_path,"/clustering_DBSCAN"))){dir.create(paste0(output_path,"/clustering_DBSCAN"))}
for(k in seq(4,7,by=.05)){
  clusters <- dbscan(DATA2@dr$tsne@cell.embeddings, eps = k ,minPts = 20)
  names(clusters$cluster) <- rownames(DATA2@meta.data)
  DATA2 <- AddMetaData(object = DATA2, metadata = clusters$cluster, col.name = paste0("dbscan.",k))
  png(filename = paste0(output_path,"/clustering_DBSCAN/tSNE_DBSCAN_k=",format(k,nsmall=2),".png"),width = 700,height = 600,res = 150)
  TSNEPlot(object = DATA2, group.by=paste0("dbscan.",k), pt.size = 1, plot.title= paste0("DBSCAN clustering with minPts=",k))
  dev.off()
}
#---------



### Plot genes of interest
#---------
rownames(DATA2@scale.data)[grep("Hba",rownames(DATA2@scale.data))]

mygenelist <- c("Tcrg-C2","Cd3e","Ncr1","CD4","Cd8a","Cd8b1","Il17a","Infg","Rorc","Cd3f","Gzmb","Gzma","Foxp3","Il10","Sell","Cd44","Cd69","Thy1","Itgam","Cd19","Cd20","Hbb-bt","Lamp1","Lamp2","Fas","Runx3","Itgal","Stat1","Ctla4","Cxcr4","Xcl1","Gata3","Il2ra","Il18")
mygenelist <- mygenelist[mygenelist %in% rownames(DATA2@scale.data)]
png(filename = paste0(output_path,"/tSNE_seleceted_genes_scaled.png"),width = 200*24,height = 200*3.5*length(unique(mygenelist))/8,res = 150)
#pdf(paste0(output_path,"/tSNE_plot_genes_per_cluster_scaled.pdf"),width = 16,height = 4*length(unique(top10$gene))/4,useDingbats = F)
FeaturePlot(object = DATA2, features.plot = unique(mygenelist), cols.use = c("grey", "blue"), reduction.use = "tsne",pt.size = 1,nCol=8)
dev.off()

#---------






### Find relevant markers for each population
#---------
DATA2@ident <- factor(NULL)
DATA2 <- SetIdent(DATA2,ident.use = DATA2@meta.data$dbscan.6.6)
DATA2 <- SubsetData(DATA2, cells.use = DATA2@cell.names[DATA2@meta.data$dbscan.6.6 > 0]) #Filter out cell with no assigned clusters
DATA2 <- BuildSNN(DATA2,reduction.type = "tsne",plot.SNN = T,k.param = 3,prune.SNN = .1)
DATA2_markers <- FindAllMarkers(object = DATA2, only.pos = T)
write.csv2(DATA2_markers,file = "marker_genes.csv",row.names = T)

png(filename = paste0(output_path,"/tSNE_after_filtering_cells_out.png"),width = 700,height = 600,res = 150)
TSNEPlot(object = DATA2,group.by="ident")
dev.off()
#---------


### Plot a heatmap with the top genes for each cluster
#---------
my_markers <- DATA2_markers[(DATA2_markers$p_val < 0.01)&(DATA2_markers$avg_logFC > 0.5),]
my_markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10

#pdf(paste0(output_path,"/Cluster_markers_heatmap.pdf"),width = 10,height = 10, useDingbats = F)
png(filename = paste0(output_path,"/Cluster_markers_heatmap.png"),width = 1300,height = 1300,res = 150)
DoHeatmap(object = DATA2, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,col.low = "grey95",col.mid = "grey85",col.high = "navy")
dev.off()
#---------




### Plot violin plots with the top genes for each cluster
#---------
png(filename = paste0(output_path,"/violinPlot_genes_per_cluster_scaled.png"),width = 200*40,height = 200*2*length(unique(top10$gene))/10,res = 150)
#pdf(paste0(output_path,"/violinPlot_genes_per_cluster_scaled.pdf"),width = 15,height = 2*length(unique(top10$gene))/4,useDingbats = F)
VlnPlot(object = DATA2, features.plot = unique(top10$gene),point.size.use = .1,nCol=10)
dev.off()

png(filename = paste0(output_path,"/tSNE_plot_genes_per_cluster_scaled.png"),width = 200*30,height = 200*3*length(unique(top10$gene))/9,res = 150)
#pdf(paste0(output_path,"/tSNE_plot_genes_per_cluster_scaled.pdf"),width = 16,height = 4*length(unique(top10$gene))/4,useDingbats = F)
FeaturePlot(object = DATA2, features.plot = unique(top10$gene), cols.use = c("grey", "blue"), reduction.use = "tsne",pt.size = 1,nCol=10)
dev.off()
#---------



### Plot a SNNgraph-connected tSNE
#---------
library(igraph)
png(filename = paste0(output_path,"/tSNE_plots/tSNE_SNN_cluster_plot.png"),width = 600,height =650,res = 100)
net <- graph.adjacency(adjmatrix = as.matrix(DATA2@snn), mode = "undirected", weighted = TRUE, diag = FALSE)
plot.igraph(x = net, layout = as.matrix(x = DATA2@dr$tsne@cell.embeddings), edge.width = E(graph = net)$weight/2, vertex.label = NA,
            edge.color = colorRampPalette(c("grey90","black"))(50)[round(E(graph = net)$weight*49+1)],
            vertex.size = 3,vertex.frame.color=hue_pal(l=50, c=80)(length(levels(DATA2@ident)))[DATA2@ident],
            vertex.color = hue_pal()(length(levels(DATA2@ident)))[DATA2@ident])
dev.off()
#---------



### Plot the percentage of metadata in each cluster and vice-versa
#---------
proportion <- as.data.frame(lapply(levels(DATA2@ident),function(x){c(unname(table(DATA2@meta.data$days_post_infection[DATA2@ident==x])[-4]))} ),row.names = c(0,3,8))
colnames(proportion) <- levels(DATA2@ident)

sa <- cbind(stack(as.data.frame(proportion/rowSums(proportion))), rep(rownames(proportion),ncol(proportion)) )
colnames(sa) <- c("prop","clusters","days_post_infection")

proportion2 <- t(as.matrix(t(proportion)/colSums(proportion)))
cl_order <- order(proportion2["0",],proportion2["3",],proportion2["8",],decreasing=T)
sa2 <- cbind(stack(as.data.frame(t(t(proportion[,cl_order])/colSums(proportion[,cl_order])) )), rep(rownames(proportion),ncol(proportion)) )
colnames(sa2) <- c("prop","clusters","days_post_infection")


png(filename = paste0(output_path,"/barplot_cluster_days_post_infection.png"),width = 5000,height = 500,res = 150)
plot_grid(ggplot(data=sa2,mapping = aes(x=clusters, y=prop, fill=days_post_infection) ) + geom_bar(stat="identity"),ncol = 8,
          ggplot(data=sa[sa$clusters==0,],mapping = aes(x=days_post_infection, y=prop, fill=days_post_infection) ) + geom_bar(stat="identity")+ ggtitle(paste0("cluster",0)),
          ggplot(data=sa[sa$clusters==1,],mapping = aes(x=days_post_infection, y=prop, fill=days_post_infection) ) + geom_bar(stat="identity")+ ggtitle(paste0("cluster",1)),
          ggplot(data=sa[sa$clusters==2,],mapping = aes(x=days_post_infection, y=prop, fill=days_post_infection) ) + geom_bar(stat="identity")+ ggtitle(paste0("cluster",2)),
          ggplot(data=sa[sa$clusters==3,],mapping = aes(x=days_post_infection, y=prop, fill=days_post_infection) ) + geom_bar(stat="identity")+ ggtitle(paste0("cluster",3)),
          ggplot(data=sa[sa$clusters==4,],mapping = aes(x=days_post_infection, y=prop, fill=days_post_infection) ) + geom_bar(stat="identity")+ ggtitle(paste0("cluster",4)),
          ggplot(data=sa[sa$clusters==5,],mapping = aes(x=days_post_infection, y=prop, fill=days_post_infection) ) + geom_bar(stat="identity")+ ggtitle(paste0("cluster",5)))
dev.off()
#---------


### Identifying relevant markers across days_post_infection for each population
#---------
marker_list <- list()
cluster_data <- list()
if(!dir.exists(paste0(output_path,"/DGE_days_post_infection"))){dir.create(paste0(output_path,"/DGE_days_post_infection"))}
for(i in unique(DATA2@ident)){
  cat(paste0("... Processing cell cluster #",i," ..."),"\n")
  temp <- SubsetData(DATA2, cells.use = DATA2@cell.names[DATA2@ident == i]) #Select cell from a cluster
  temp@ident <- factor(NULL)
  temp <- SetIdent(temp,ident.use = temp@meta.data$days_post_infection)
  temp_markers <- FindAllMarkers(object = temp, only.pos = T)
  temp_markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top_temp
  temp_markers <- temp_markers[(temp_markers$p_val < 0.01)&(temp_markers$avg_logFC > 0.3),]
  temp_markers <- temp_markers[order(temp_markers$p_val),]
  marker_list[[i]] <- temp_markers
  cluster_data[[i]] <- as.matrix(temp@data)[temp_markers$gene,]
  write.csv2(temp_markers,paste0(output_path,"/DGE_days_post_infection/markers_embryonic.age_for_cluster",i,".csv"),row.names = T)
  if(nrow(temp_markers) > 20){top_temp <- temp_markers[1:20,]}else{top_temp <- temp_markers}
  
  write.csv2(temp_markers,paste0(output_path,"/DGE_days_post_infection/markers_embryonic.age_for_cluster",i,".csv"),row.names = T)
  png(filename = paste0(output_path,"/DGE_days_post_infection/vioplot_for_cluster",i,".png"),width = 200*5,height = 200*1.5*length(unique(top_temp$gene))/4,res = 150)
  print(VlnPlot(object = temp, features.plot = unique(top_temp$gene), point.size.use = .1)) #it does not work if you don't have the print command in front of it!
  dev.off()
}
#---------


for(i in unique(DATA2@ident)){
  cat(paste0("... Processing cell cluster #",i," ..."),"\n")
  temp <- SubsetData(DATA2, cells.use = DATA2@cell.names[DATA2@ident == i]) #Select cell from a cluster
  temp@ident <- factor(NULL)
  temp <- SetIdent(temp,ident.use = temp@meta.data$days_post_infection)
  top_temp <- temp_markers[(temp_markers$p_val < 0.001)&(temp_markers$avg_logFC > 0.5),]
  png(filename = paste0(output_path,"/DGE_days_post_infection/vioplot_for_cluster",i,".png"),width = 200*5,height = 200*1.5*length(unique(top_temp$gene))/4,res = 150)
  print(VlnPlot(object = temp, features.plot = unique(top_temp$gene), point.size.use = .1)) #it does not work if you don't have the print command in front of it!
  dev.off()
}










### Plotting a heatmap for differentially expressed genes among time points
for(i in unique(DATA2@ident)){
  annotation_col = data.frame( Embryonic_age = DATA2@meta.data[colnames(cluster_data[[i]]),"days_post_infection"],row.names = colnames(cluster_data[[i]]))
  ann_colors = list(  Embryonic_age = hue_pal()(length(unique(DATA2@meta.data[,"days_post_infection"]))) )
  names(ann_colors$Embryonic_age) <- levels(DATA2@meta.data[,"days_post_infection"])
  if(nrow(cluster_data[[i]]) > 30){sel <- 30} else {sel <- nrow(cluster_data[[i]])}
  pheatmap(cluster_data[[i]][1:sel,order(DATA2@meta.data[colnames(cluster_data[[i]]),"days_post_infection"])],scale = "row",color = colorRampPalette(c("navy","navy","white","red","red"))(50),
           border_color = NA,clustering_method = "ward.D2",annotation_col = annotation_col, annotation_colors = ann_colors,
           cluster_cols = F,gaps_col = cumsum(table(DATA2@meta.data[colnames(cluster_data[[i]]),"days_post_infection"]))[-3],
           filename=paste0(output_path,"/DGE_days_post_infection/heatmap_for_cluster",i,".png"))
}
#---------



### Saving the Seurat object
#---------
saveRDS(DATA2, file = paste0(output_path,"/T_cell_Seurat_object.rds") )
#---------



### Run Monocle 2 using the tSNE space
#---------
source("http://bioconductor.org/biocLite.R")
biocLite("monocle")
library(monocle)
if(!dir.exists(paste0(output_path,"/monocle2"))){dir.create(paste0(output_path,"/monocle2"))}


newtemp <- newCellDataSet(as.matrix(DATA2@scale.data), phenoData = new("AnnotatedDataFrame", data = cbind(DATA2@meta.data,ident=factor(DATA2@ident) ) ) ,featureData = NULL)
newtemp <- estimateSizeFactors(newtemp)

#diff_test_res <- differentialGeneTest(newtemp, fullModelFormulaStr = "~ident",relative_expr = F)
#ordering_genes <- row.names(subset(diff_test_res[order(diff_test_res$qval),][diff_test_res$status=="OK",], qval < 0.01))
#newtemp <-  setOrderingFilter(newtemp,  ordering_genes = ordering_genes)

my_markers %>% group_by(cluster) %>% top_n(20, avg_logFC) -> top20
newtemp <-  setOrderingFilter(newtemp,  ordering_genes = top20$gene)

#ordering_genes2 <- rownames(DATA2_markers[(DATA2_markers$p_val_adj<0.01) & (DATA2_markers$avg_logFC>0.5),])
#newtemp <-  setOrderingFilter(newtemp,  ordering_genes = ordering_genes2)

newtemp <- reduceDimension(newtemp,reduction_method="ICA")
newtemp@reducedDimS <- t(DATA2@dr$tsne@cell.embeddings)

newtemp <- orderCells(newtemp)
plot_cell_trajectory(newtemp, color_by = "State")
plot_cell_trajectory(newtemp, color_by = "embrionic.age")
plot_cell_trajectory(newtemp, color_by = "ident")
newtemp <- orderCells(newtemp,root_state = 5,num_paths = 4)

png(filename = paste0(output_path,"/monocle2/monocle_trajectories.png"),width = 3.5*700,height = 2*800,res = 150)
plot_grid(plot_cell_trajectory(newtemp, color_by = "State"),
          plot_cell_trajectory(newtemp, color_by = "embrionic.age"),
          plot_cell_trajectory(newtemp, color_by = "Pseudotime"),
          plot_cell_trajectory(newtemp, color_by = "ident"),
          plot_cell_trajectory(newtemp, color_by = "res.0.3"), ncol = 3)
dev.off()


png(filename = paste0(output_path,"/monocle2/monocle_per_trajectories.png"),width = 5*700,height = 3*800,res = 150)
plot_grid(plot_cell_trajectory(newtemp, color_by = "State") + facet_wrap(~State, nrow = 1),
          plot_cell_trajectory(newtemp, color_by = "embrionic.age") + facet_wrap(~embrionic.age, nrow = 1),
          plot_cell_trajectory(newtemp, color_by = "ident") + facet_wrap(~ident, nrow = 1),
          ncol=1)
dev.off()

DATA2 <- AddMetaData(object = DATA2, metadata = as.matrix(pData(newtemp))[,"State"], col.name = "Monocle.State")
DATA2 <- AddMetaData(object = DATA2, metadata = as.matrix(pData(newtemp))[,"Pseudotime"], col.name = "Monocle.Pseudotime")
DATA2@meta.data$Monocle.Pseudotime <- as.numeric(as.matrix(pData(newtemp))[,"Pseudotime"])
png(filename = paste0(output_path,"/monocle2/tSNE_monocle_trajectories.png"),width = 3*700,height = 700,res = 150)
plot_grid(TSNEPlot(object = DATA2, do.return = TRUE, group.by = c("ident"), no.legend = TRUE, do.label = TRUE) + ggtitle(paste0("cell ident")),
          #FeaturePlot(object = DATA2, features.plot = "Monocle.Pseudotime",cols.use = col_scale, reduction.use = "tsne") + ggtitle(paste0("Monocle2 pseudotime")),
          TSNEPlot(object = DATA2, do.return = TRUE, group.by = c("Monocle.State"), no.legend = TRUE, do.label = TRUE) + ggtitle(paste0("Monocle2 cell state")),
          TSNEPlot(object = DATA2, do.return = TRUE, group.by = c("res.0.3"), no.legend = TRUE, do.label = TRUE) + ggtitle(paste0("res0.3")),
          ncol = 3)
dev.off()


png(filename = paste0(output_path,"/monocle2/tSNE_monocle_pseudotime.png"),width = 700,height = 700,res = 150)
FeaturePlot(object = DATA2, features.plot = "Monocle.Pseudotime",cols.use = rev(col_scale) )
dev.off()
#---------



png(filename = paste0(output_path,"/monocle2/Pseudotime_per_group.png"),width = 2100,height = 700,res = 150)
plot_grid(VlnPlot(DATA2,features.plot = "Monocle.Pseudotime",do.sort = T),
          VlnPlot(DATA2,features.plot = "Monocle.Pseudotime",group.by = "embrionic.age"),
          VlnPlot(DATA2,features.plot = "Monocle.Pseudotime",group.by = "collection_date"),ncol = 3)
dev.off()
#---------





mylist <- c("Bmp4","Nog","Alpl","Ngfr","Wnt5a","Lef1","Ptc1","Gli1","Pfgfra","Dkk1","Dkk2","Dkk4","Bmp2","Bmp4","Bmp7","Col6a1","Col6a5","Cd36","Wif1","Apcdd1","Hspb3","Sox9")
mylist[!(mylist %in% rownames(DATA2@data))]
mylist <- mylist[mylist %in% rownames(DATA2@data)]
png(filename = paste0(output_path,"/tSNE_genes_of_interest.png"),width = 2100,height = 2600,res = 150)
FeaturePlot(object = DATA2, features.plot = mylist,cols.use = col_scale )
dev.off()






### Run cellAlign
#---------
install_github("shenorrLab/cellAlign")


#---------








#if(!dir.exists(paste0(output_path,"/clustering_SNN"))){dir.create(paste0(output_path,"/clustering_SNN"))}
#p <- list()
#for(k in seq(.05,4,by=.05)){
#  if(k!=.05){
#     DATA2 <- FindClusters(object = DATA2, reduction.type = "pca", dims.use = 1:10, resolution = k, print.output = F)
#   } else { DATA2@ident <- factor(NULL) ; DATA2 <- FindClusters(object = DATA2, reduction.type = "pca", dims.use = 1:10, resolution = k, print.output = F, save.SNN = TRUE,force.recalc = T)}
#   #png(filename = paste0(output_path,"/clustering_SNN/tSNE_SNN_res=",k,".png"),width = 700,height = 600,res = 150)
#   p[[paste0(k)]] <- TSNEPlot(object = DATA2, group.by=paste0("res.",k), pt.size = 1, plot.title= paste0("SNN clustering with res=",k))
#   #dev.off()
# }
# 
# png(filename = paste0(output_path,"/clustering_SNN/tSNE_SNN.png"),width = 4000,height = 3000,res = 100)
# do.call(plot_grid,p)
# dev.off()


# p <- list(all=TSNEPlot(object = DATA2,group.by="embrionic.age",pt.size = 1,colors.use=rep("grey",30) ),
#           per_time=TSNEPlot(object = DATA2,group.by="embrionic.age",pt.size = 1) )
# for(j in unique(DATA2@meta.data$embrionic.age)){
#   a <- SubsetData(DATA2, cells.use = DATA2@cell.names[DATA2@meta.data$embrionic.age == j])
#   p[[j]] <- TSNEPlot(object = a,group.by="embrionic.age",pt.size = 1,colors.use = hue_pal()(length(unique(DATA2@meta.data$embrionic.age)))[unique(DATA2@meta.data$embrionic.age)[unique(DATA2@meta.data$embrionic.age)==j]]) }
# png(filename = paste0(output_path,"/tSNE_plots/tSNE_per_embryonic_age.png"),width = 1800,height = 1000,res = 150)
# do.call(plot_grid,p)
# rm(p)
# dev.off()





mito.genes <- grep(pattern = "^mt-", x = rownames(x = DATA2@data), value = TRUE)
percent.mito <- Matrix::colSums(DATA2@raw.data[mito.genes, ]) / Matrix::colSums(DATA2@raw.data)
DATA2 <- AddMetaData(object = DATA2, metadata = percent.mito, col.name = "percent.mito")

Ubq.genes <- grep(pattern = "^Ub.[123456789]", x = rownames(x = DATA2@data), value = TRUE)
percent.Ubq <- Matrix::colSums(DATA2@raw.data[Ubq.genes, ]) / Matrix::colSums(DATA2@raw.data)
DATA2 <- AddMetaData(object = DATA2, metadata = percent.Ubq, col.name = "percent.Ubq")

Rps.genes <- grep(pattern = "^Rps[123456789]", x = rownames(x = DATA2@data), value = TRUE)
percent.Rps <- Matrix::colSums(DATA2@raw.data[Rps.genes, ]) / Matrix::colSums(DATA2@raw.data)
DATA2 <- AddMetaData(object = DATA2, metadata = percent.Rps, col.name = "percent.Rps")

Rpl.genes <- grep(pattern = "^Rpl[123456789]", x = rownames(x = DATA2@data), value = TRUE)
percent.Rpl <- Matrix::colSums(DATA2@raw.data[Rpl.genes, ]) / Matrix::colSums(DATA2@raw.data)
DATA2 <- AddMetaData(object = DATA2, metadata = percent.Rpl, col.name = "percent.Rpl")

Zf.genes <- grep(pattern = "^Zfp", x = rownames(x = DATA2@data), value = TRUE)
percent.Zf <- Matrix::colSums(DATA2@raw.data[Zf.genes, ]) / Matrix::colSums(DATA2@raw.data)
DATA2 <- AddMetaData(object = DATA2, metadata = percent.Zf, col.name = "percent.Zf")

FeaturePlot(object = DATA2, features.plot = c("percent.mito","percent.Rpl","percent.Rps","percent.Ubq","percent.Zf"),cols.use = col_scale )

plot(percent.Rps,percent.Rpl)


Gene.groups <- substring(rownames(x = DATA2@data),1,3)
temp <- rowsum(as.matrix(DATA2@raw.data),Gene.groups) / colSums(DATA2@raw.data)
perc <- sort(rowMeans(temp),decreasing = T)

mypar(2,1)
boxplot(100*t(temp[rownames(temp)%in%names(perc)[1:40],])[,names(perc)[1:40]],outline=F,las=2,ylab="% reads")
barplot(perc[1:40]*100,las=2,xaxs="i",ylab="% reads")





