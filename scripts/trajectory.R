




DATA <- readRDS("~/Downloads/N_Lycke_1801/analysis/2_clustering_batchDPI/Seurat_object.rds")
DATA2 <- RunUMAP(DATA2,verbose=T,max.dim = 10,min_dist = 0.00001,n_neighbors = 80,dims.use = 1:20)
DimPlot(DATA2,reduction.use = "umap",group.by = "res.0.5",dim.1 = 1,dim.2 = 2)
FeaturePlot(object = DATA2, features.plot = "Krt14", cols.use = c("grey","blue"), reduction.use = "umap",pt.size = 1,dim.1 = 1,dim.2 = 2)
FeaturePlot(object = DATA2, features.plot = "Krt17", cols.use = c("grey","blue"), reduction.use = "umap",pt.size = 1,dim.1 = 1,dim.2 = 2)
FeaturePlot(object = DATA2, features.plot = "Dkk4", cols.use = c("grey","blue"), reduction.use = "umap",pt.size = 1,dim.1 = 1,dim.2 = 2)

#Path using UMAP
#----------------
newtemp <- newCellDataSet(cellData = as.matrix(DATA2@raw.data)[,], 
                          phenoData = new("AnnotatedDataFrame", data = cbind(DATA2@meta.data[,] ) ) ,
                          featureData = NULL)
newtemp <- estimateSizeFactors(newtemp)

#newtemp <- reduceDimension(newtemp,reduction_method="ICA",max_components = 2)
newtemp@reducedDimS <- t(DATA2@dr$umap@cell.embeddings[,])
newtemp@reducedDimW <- 0; newtemp@reducedDimA <- 0; newtemp@reducedDimK <- 0

newtemp <- orderCells(newtemp,num_paths = 4)
newtemp@phenoData@data$State <- newtemp@phenoData@data$res.0.4
newtemp <- orderCells(newtemp,root_state = "6",num_paths = 4)

plot_cell_trajectory(newtemp, color_by = "State")
plot_cell_trajectory(newtemp, color_by = "res.0.4")
plot_cell_trajectory(newtemp, color_by = "Pseudotime")
#----------------






















library(scales)

DATA <- RunUMAP(DATA,verbose=T,max.dim = 2,min_dist = 0.00001,n_neighbors = 50,dims.use = 1:50)
DATA2 <- RunUMAP(DATA,verbose=T,max.dim = 10,min_dist = 0.00001,n_neighbors = 50,dims.use = 1:50)

d <- dist(DATA2@dr$umap@cell.embeddings,method = "euclidean")
hc <- hclust(d,method = "ward.D2")
n_c <- round(sqrt(nrow(DATA@dr$umap@cell.embeddings)))
cl <- cutree(hc,k = n_c)
plot(as.data.frame(DATA@dr$umap@cell.embeddings)[names(cl),],col=hue_pal()(n_c)[factor(cl)],cex=.4,pch=16)


#d <- as.dist( -(1-cor(t(DATA2@dr$umap@cell.embeddings)))/2  )
#hc <- hclust(d,method = "ward.D2")
#n_c <- round(sqrt(nrow(DATA@dr$umap@cell.embeddings))/2)
#cl <- cutree(hc,k = n_c)
#plot(as.data.frame(DATA@dr$umap@cell.embeddings)[names(cl),],col=hue_pal()(n_c)[factor(cl)],cex=.4,pch=16)


temp <- rowsum(DATA@dr$umap@cell.embeddings,group = cl)
temp <- temp / c(table(cl))
points(temp,col="#FFFFFF99",cex=2,pch=16)
text(temp,labels = rownames(temp),cex=.8)


temp <- rowsum(DATA2@dr$umap@cell.embeddings,group = factor(cl))
temp <- temp / c(table(cl))
d_c <- dist(temp,method = "euclidean")
hc_c <- hclust(d_c,method = "ward.D2")
plot(as.dendrogram(hc_c),las=1)
points( 1:length(hc_c$order),rep(0,length(hc_c$order)),col= hue_pal()(n_c)[factor(hc_c$order)] ,pch=16,cex=2)










plot(rnorm(n = length(y),sd = 0.1),y,xlim=c(-1,5),cex=.1)
points(rnorm(n = length(y2),sd = 0.1)+1,y2,cex=.1)
points(rnorm(n = length(y2),sd = 0.1)+2,DATA@data[,3],cex=.1)
points(rnorm(n = length(y2),sd = 0.1)+3,DATA@data[,4],cex=.1)
points(rnorm(n = length(y2),sd = 0.1)+4,DATA@data[,5],cex=.1)


library(monocle)
library(Seurat)
opt <- list()
opt$clustering_use <- "res.0.4"


DATA <- readRDS("~/Downloads/N_Lycke_1801/analysis/2_clustering_batchDPI/Seurat_object.rds")
DATA <- RunUMAP(DATA,verbose=T,max.dim = 2,min_dist = 0.00001,n_neighbors = 30,dims.use = 1:20)
DimPlot(DATA,reduction.use = "umap",group.by = "res.0.5")
DimPlot(DATA,reduction.use = "umap",group.by = "days_post_infection")



temp2 <- FindClusters(temp,reduction.type = "umap",force.recalc = T,resolution = .5,algorithm = 3)
DimPlot(temp2,reduction.use = "umap",group.by = "res.0.5")



temp <- RunUMAP(DATA,verbose=T,max.dim = 3,min_dist = 0.00001,n_neighbors = 50,dims.use = 1:20)
DimPlot(temp,reduction.use = "umap",group.by = "res.0.5",dim.1 = 1,dim.2 = 2)
DimPlot(temp,reduction.use = "umap",group.by = "res.0.5",dim.1 = 1,dim.2 = 3)
DimPlot(temp,reduction.use = "umap",group.by = "days_post_infection",dim.1 = 1,dim.2 = 3)


### Run Monocle 2 using the UMAP space
#---------
sel <- colnames(temp@raw.data)[sample(x = 1:(dim(temp@raw.data)[2]),size = 2000)]
newtemp <- newCellDataSet(cellData = as.matrix(temp@raw.data)[,sel], 
                          phenoData = new("AnnotatedDataFrame", data = cbind(temp@meta.data[sel,] ) ) ,
                          featureData = NULL)
newtemp <- estimateSizeFactors(newtemp)



if(opt$genes_use == "monocle"){
  diff_test_res <- differentialGeneTest(newtemp, fullModelFormulaStr = paste0("~",opt$clustering_use),relative_expr = F)
  ordering_genes <- row.names(subset(diff_test_res[order(diff_test_res$qval),][diff_test_res$status=="OK",], qval < 0.01))
  newtemp <-  setOrderingFilter(newtemp,  ordering_genes = ordering_genes)
} else if (opt$genes_use == "var.genes"){
  newtemp <-  setOrderingFilter(newtemp,  ordering_genes = temp@var.genes)
} else if(opt$genes_use == "marker.genes"){
  markers <- read.csv2("~/Downloads/N_Lycke_1801/analysis/2_clustering_batchDPI/diff_expr/Cluster_marker_genes.csv")
}


rownames(diff_test_res)[ diff_test_res$use_for_ordering==TRUE ]

#ordering_genes2 <- rownames(DATA2_markers[(DATA2_markers$p_val_adj<0.01) & (DATA2_markers$avg_logFC>0.5),])
#newtemp <-  setOrderingFilter(newtemp,  ordering_genes = ordering_genes2)

dispersionTable(newtemp)





DATA <- RunICA(DATA)
DimPlot(DATA,reduction.use = "ica",group.by = "res.0.4",dim.1 = 1,dim.2 = 2)
DimPlot(DATA,reduction.use = "ica",group.by = "res.0.4",dim.1 = 3,dim.2 = 4)
DimPlot(DATA,reduction.use = "ica",group.by = "res.0.4",dim.1 = 5,dim.2 = 6)


DATA <- RunDiffusion(DATA, dims.use = 1:10, max.dim = 10)
DimPlot(DATA,reduction.use = "dm",group.by = "res.0.4",dim.1 = 1,dim.2 = 2)
DimPlot(DATA,reduction.use = "dm",group.by = "res.0.4",dim.1 = 3,dim.2 = 4)
DimPlot(DATA,reduction.use = "dm",group.by = "res.0.4",dim.1 = 5,dim.2 = 6)
DimPlot(DATA,reduction.use = "dm",group.by = "res.0.4",dim.1 = 7,dim.2 = 8)
DimPlot(DATA,reduction.use = "dm",group.by = "res.0.4",dim.1 = 9,dim.2 = 10)








#Path using DDRTree
#----------------
newtemp3 <- reduceDimension(newtemp,reduction_method="DDRTree")
newtemp3 <- orderCells(newtemp3)
plot_cell_trajectory(newtemp3, color_by = "res.0.4")
plot_cell_trajectory(newtemp3, color_by = "Pseudotime")
plot_cell_trajectory(newtemp3, color_by = "State")
newtemp@phenoData@data$State <- newtemp@phenoData@data$res.0.4

newtemp <- orderCells(newtemp,root_state = "14")
plot_cell_trajectory(newtemp, color_by = "res.0.4")
plot_cell_trajectory(newtemp, color_by = "Pseudotime")
#----------------





#Path using ICA
#----------------
newtemp <- reduceDimension(newtemp,reduction_method="ICA",max_components = 2)
newtemp <- orderCells(newtemp,num_paths = 4)
newtemp@phenoData@data$State <- newtemp@phenoData@data$res.0.4
newtemp <- orderCells(newtemp,root_state = "3",num_paths = 4)
plot_cell_trajectory(newtemp, color_by = "res.0.4")
plot_cell_trajectory(newtemp, color_by = "Pseudotime")
plot_cell_trajectory(newtemp, color_by = "State")
#----------------





#Path using Difusion maps
#----------------
newtemp2 <- reduceDimension(newtemp,reduction_method="ICA",max_components = 10)
newtemp2@reducedDimS <- t(temp@dr$dm@cell.embeddings[sel,1:10])
newtemp2 <- orderCells(newtemp2,num_paths = 1)
newtemp2@phenoData@data$State <- newtemp2@phenoData@data$res.0.4
newtemp2 <- orderCells(newtemp2,root_state = "3",num_paths = 4)

plot_cell_trajectory(newtemp2, color_by = "State")
plot_cell_trajectory(newtemp2, color_by = "days_post_infection")
plot_cell_trajectory(newtemp2, color_by = "Pseudotime")
#----------------






#Path using UMAP
#----------------
newtemp2 <- reduceDimension(newtemp,reduction_method="ICA",max_components = 2)
newtemp2@reducedDimS <- t(DATA@dr$umap@cell.embeddings[sel,])

newtemp2 <- orderCells(newtemp2,num_paths = 4)
newtemp2@phenoData@data$State <- newtemp2@phenoData@data$res.0.4
newtemp2 <- orderCells(newtemp2,root_state = "3",num_paths = 4)

plot_cell_trajectory(newtemp2, color_by = "State")
plot_cell_trajectory(newtemp2, color_by = "res.0.4")
plot_cell_trajectory(newtemp2, color_by = "Pseudotime")
#----------------





temp2 <- SubsetData(DATA, cells.use = sel) #Select cell from a cluster
temp2@meta.data[sel,"State"] <- newtemp2@phenoData@data[sel,"State"]
temp2@meta.data[sel,"Pseudotime"] <- newtemp2@phenoData@data[sel,"Pseudotime"]
DimPlot(temp2,reduction.use = "umap",group.by = "State")
DimPlot(temp2,reduction.use = "umap",group.by = "res.0.4")
FeaturePlot(object = temp2, features.plot = "Pseudotime", cols.use = colorRampPalette( c("grey10", "grey","blue") )(90), reduction.use = "umap",pt.size = 1)
FeaturePlot(object = temp2, features.plot = "Pseudotime", cols.use = colorRampPalette( c("grey10", "grey","blue") )(90), reduction.use = "tsne",pt.size = 1)









diff_test_res_pseudotime <- differentialGeneTest(newtemp, fullModelFormulaStr = "~sm.ns(Pseudotime)")
top_pseudotime <-rownames(diff_test_res_pseudotime[order(diff_test_res_pseudotime$qval,decreasing = F),])[1:500]
plot_pseudotime_heatmap(newtemp[top_pseudotime,], num_clusters = 7, cores = 4, show_rownames = T)





newtemp





my_genes <- c("Il17a","Sell","Nkg7","Gzmk","Cd44","Ifng")
newtemp_subset <- newtemp[my_genes,]
plot_genes_branched_pseudotime(newtemp, branch_point = 1, color_by = "Pseudotime", ncol = 1)
plot_genes_branched_pseudotime(newtemp_subset, branch_point = 1, color_by = "res.0.4", ncol = 1)
plot_cell_trajectory(newtemp, color_by = "Pseudotime")
plot_genes_in_pseudotime(newtemp_subset, color_by = "res.0.4")


plot_cell_trajectory(newtemp, markers = c("Cd44", "Sell"), use_color_gradient = TRUE)

diff_test_res <- differentialGeneTest(newtemp_subset, fullModelFormulaStr = "~sm.ns(Pseudotime)")


temp2 <- SubsetData(DATA, cells.use = sel) #Select cell from a cluster
temp2@meta.data[sel,"State"] <- newtemp@phenoData@data[sel,"State"]
temp2@meta.data[sel,"Pseudotime"] <- newtemp@phenoData@data[sel,"Pseudotime"]
DimPlot(temp2,reduction.use = "umap",group.by = "State")
DimPlot(temp2,reduction.use = "umap",group.by = "res.0.4")
FeaturePlot(object = temp2, features.plot = "Pseudotime", cols.use = colorRampPalette( c("grey10", "grey","blue") )(90), reduction.use = "umap",pt.size = 1)
FeaturePlot(object = temp2, features.plot = "Pseudotime", cols.use = colorRampPalette( c("grey10", "grey","blue") )(90), reduction.use = "tsne",pt.size = 1)


FeaturePlot(object = temp2, features.plot = "Itgal", cols.use = c("grey", "navy"), reduction.use = "umap",pt.size = 1)




#Path using UMAP
#----------------
newtemp <- reduceDimension(newtemp,reduction_method="ICA",max_components = 10)
dim(newtemp@reducedDimS)
#newtemp@reducedDimS <- t(DATA@dr$dm@cell.embeddings[sel,])


newtemp <- orderCells(newtemp,num_paths = 1)
plot_cell_trajectory(newtemp, color_by = "res.0.4")
plot_cell_trajectory(newtemp, color_by = "Pseudotime")
plot_cell_trajectory(newtemp, color_by = "State")

newtemp@phenoData@data$State <- newtemp@phenoData@data$res.0.4


newtemp <- orderCells(newtemp,root_state = "3",num_paths = 4)
plot_cell_trajectory(newtemp, color_by = "res.0.4")
plot_cell_trajectory(newtemp, color_by = "Pseudotime")


plot_cell_trajectory(newtemp, color_by = "Pseudotime")
rownames(newtemp$Pseudotime)

plot_cell_trajectory(newtemp, color_by = "State")


temp2 <- SubsetData(DATA, cells.use = DATA@cell.names[sel]) #Select cell from a cluster
temp2@meta.data$Pseudotime <- newtemp$Pseudotime
DimPlot(temp2,reduction.use = "umap",group.by = "res.0.4")
FeaturePlot(object = temp2, features.plot = "Pseudotime", cols.use = c("grey", "navy"), reduction.use = "umap",pt.size = 1)
FeaturePlot(object = temp2, features.plot = "Pseudotime", cols.use = c("grey", "navy"), reduction.use = "tsne",pt.size = 1)
#----------------











#Path using UMAP
#----------------
newtemp2 <- reduceDimension(newtemp,reduction_method="ICA",max_components = 5)
newtemp2@reducedDimS <- t(temp@dr$umap@cell.embeddings[sel,1:5])
newtemp2 <- orderCells(newtemp2,num_paths = 1)
newtemp2@phenoData@data$State <- newtemp2@phenoData@data$res.0.4
newtemp2 <- orderCells(newtemp2,root_state = "3",num_paths = 4)

plot_cell_trajectory(newtemp2, color_by = "State")
plot_cell_trajectory(newtemp2, color_by = "days_post_infection")
plot_cell_trajectory(newtemp2, color_by = "Pseudotime")
#----------------









newtemp <- orderCells(newtemp,root_state = 5,num_paths = 4,)

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




