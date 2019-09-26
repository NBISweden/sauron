#!/usr/bin/env Rscript
################################
### DEFINE SCRIPT PARAMETERS ###
################################
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')}
library(optparse)
cat("\nRunning TRAJECTORY ANALYSIS with the following parameters ...\n")
option_list = list(
make_option(c("-i", "--Seurat_object_path"),    type = "character",   metavar="character",   default='none',  help="Path to the Seurat object FILE."),
make_option(c("-m", "--metadata_use"),          type = "character",   metavar="character",   default='none',  help="Column names of the metadata table to plot in the trajectory map."),
make_option(c("-d", "--reduction_use"),         type = "character",   metavar="character",   default='umap10',  help="Dimensionality reduction method to base the trajectories on. It could be a pre-computed slot within your Seurat Object (in case: pca, umap, umap10, tsne) or it will compute additional others (ica, dm). Trajectories are defined on top of UMAP ambedding with 10 dimensions by default."),
make_option(c("-r", "--reduction_visualize"),   type = "character",   metavar="character",   default='umap',  help="Dimensionality reduction method to visualize trajectories. It could be a pre-computed slot within your Seurat Object (in case: pca, umap, umap10, tsne) or it will compute additional others (ica, dm). Trajectories are defined on top of UMAP ambedding with 10 dimensions by default."),
make_option(c("-e", "--method_use"),            type = "character",   metavar="character",   default='DDRTree',  help="Method used for trajectory estimation, MST or DDRTree. DDRTree is default."),
make_option(c("-t", "--no_traj_components"),    type = "character",   metavar="character",   default='3',  help="Number trajectory components to use"),
make_option(c("-p", "--no_of_paths"),           type = "character",   metavar="character",   default='none',  help="Number specifying how many paths to allow."),
make_option(c("-n", "--cluster_use"),           type = "character",   metavar="character",   default='all',    help="The cluster of cells to be used for analysis. Should be defined as the clustering name followed by the cluster names to be used, comma-separated. E.g.: 'louvain_0.2,1,2,3,5,6'."),
make_option(c("-s", "--start_cluster"),         type = "character",   metavar="character",   default='none',  help="Cluster from which the trajectories will start from."),
make_option(c("-z", "--diff_testing"),          type = "character",   metavar="character",   default='none',  help="Whether to test for diffential expression across branches"),
make_option(c("-a", "--assay"),                 type = "character",   metavar="character",   default='RNA',  help="Assay to use for trajectory differential expression. Default is 'RNA'"),
make_option(c("-o", "--output_path"),           type = "character",   metavar="character",   default='none',  help="Output DIRECTORY.")
)
opt = parse_args(OptionParser(option_list=option_list))
print(t(t(unlist(opt))))

if(!dir.exists(opt$output_path)){dir.create(opt$output_path,recursive = T)}
setwd(opt$output_path)
#---------


# opt <- list(
#   Seurat_object_path = "/Users/paulo.barenco/Desktop/NBIS/Projects/M_Kasper_1709/trajectory_analysis_epithelial_190917/analysis/integrated_1000_k10/clustering/seurat_object.rds",
#   metadata_use       = "louvain_0.9"  ,
#   reduction_use      = "umap10"  ,
#   reduction_visualize= "umap"    ,
#   method_use         = "ddrtree" ,
#   no_traj_components = "5"   ,
#   no_of_paths        = "6"  ,
#   cluster_use        = "none",
#   start_cluster      = "none",
#   diff_testing       = "no"  ,
#   assay              = "RNA" ,
#   output_path       =  "/Users/paulo.barenco/Desktop/NBIS/Projects/M_Kasper_1709/trajectory_analysis_epithelial_190917/analysis/integrated_1000_k10/clustering/traj_ddrtree"
# 
# 
# )


##############################
### LOAD/INSTALL LIBRARIES ###
##############################
cat("\nLoading/installing libraries ...\n")
initial.options <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(sub("--file=","",initial.options[grep("--file=",initial.options)]))
source( paste0(script_path,"/inst_packages.R") )
pkgs <- c("Seurat","scales","monocle","plot3D","scales","parallel")
inst_packages(pkgs)
#---------



#############################
### LOAD Seurat.v3 OBJECT ###
#############################
cat("\n### LOADING Seurat.v3 OBJECT ###\n")
DATA <- readRDS(opt$Seurat_object_path)
#---------



###################################
### SELECT CELLS FROM A CLUSTER ###
###################################
cat("\n### SELECTING CELLS FROM A CLUSTER ###")
if (length(unlist(strsplit(opt$cluster_use,","))) >= 2 ){
    clustering_use <- as.character(unlist(strsplit(opt$cluster_use,",")))[1]
    clusters_to_select <- as.character(unlist(strsplit(opt$cluster_use,",")))[-1]
    if(!(clustering_use %in% colnames(DATA@meta.data))){
        cat("\nThe Clustering level was not found in the Metadata...\n")
        cat("\nAll cells will be used ...\n")
        cells_use <- rownames(DATA@meta.data)
    } else if(sum(clusters_to_select %in% unique(DATA@meta.data[,clustering_use])) == 0){
        cat("\nThe Cluster specifed was not found ...\n")
        cat("\nAll cells will be used ...\n")
        cells_use <- rownames(DATA@meta.data)
    } else {
        cells_use <- rownames(DATA@meta.data)[factor(DATA@meta.data[,clustering_use]) %in% clusters_to_select]}   #Filter out cells with no assigned clusters
    DATA <- SubsetData(DATA,assay = opt$assay,cells = cells_use)
} else {
    cat("\nThe name of the cluster or the cluster name were not found in your data.\n All cells will be used ...\n")
    cells_use <- rownames(DATA@meta.data) }
# sel <- rowSums(as.matrix(DATA@assays$RNA@counts) >= 1) >= 1
# DATA <- CreateSeuratObject(as.matrix(DATA@assays$RNA@counts[sel,cells_use]), meta.data = DATA@meta.data[cells_use,])
#---------



#########################
### DEFINE PARAMETERS ###
#########################
reduction_use <- unlist(strsplit( opt$reduction_use , split = ","))
reduction_use <- reduction_use[reduction_use %in% names(DATA@reductions)]
cat("\nThe folowing reductions were found in your dataset:\n")
print(reduction_use)
red_vis <- unlist(strsplit( opt$reduction_visualize , split = ","))
#---------



for( reduc in reduction_use){
    cat('\n\n PROCESSING reduction: ',reduc,'  ...')
    if(!dir.exists(paste0(opt$output_path,"/",reduc))){dir.create(paste0(opt$output_path,"/",reduc),recursive = T)}
    
    if(file.exists(paste0(opt$output_path,"/",reduc,"/Monocle_trajectory_object.rds"))){
        tempUMAP3 <- readRDS(paste0(opt$output_path,"/",reduc,"/Monocle_trajectory_object.rds"))
    } else {
        
        ########################
        ### BUILD SCE OBJECT ###
        ########################
        cat("\nConstructing monocle object ...\n")
        #For the whole dataset as monocle object
        temp <- newCellDataSet(cellData = as.matrix(DATA@assays[[opt$assay]]@counts[rownames(DATA@assays[[opt$assay]]@counts),]),
        phenoData = new("AnnotatedDataFrame", data = cbind(DATA@meta.data ) ) ,
        featureData = new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(DATA@assays[[opt$assay]]@counts),row.names = rownames(DATA@assays[[opt$assay]]@counts))) )
        temp <- estimateSizeFactors(temp)
        
        #For the whole data
        if( casefold(opt$method_use) == "ddrtree"){
            
            #temp2 <- newCellDataSet(cellData = t(DATA@reductions[[reduc]]@cell.embeddings),
            #phenoData = new("AnnotatedDataFrame", data = cbind(DATA@meta.data ) ) )
            
            temp2 <- newCellDataSet(cellData = DATA@assays$mnn@scale.data,
                                    phenoData = new("AnnotatedDataFrame", data = cbind(DATA@meta.data ) ) )
            
            temp2 <- estimateSizeFactors(temp2)
            tempUMAP2 <- reduceDimension(temp2, reduction_method = "DDRTree",norm_method = "none",max_components = as.numeric(opt$no_traj_components))
        } else if(casefold(opt$method_use) == "mst"){
            tempUMAP2 <- reduceDimension( temp, reduction_method = "ICA",norm_method = "none" )
            tempUMAP2@reducedDimS <- t(DATA@reductions[[reduc]]@cell.embeddings)
        }
        #---------
        
        
        
        
        ###################
        ### ORDER CELLS ###
        ###################
        cat("\nOrdering cells ...\n")
        tempUMAP2 <- orderCells(tempUMAP2)
        
        orig.state <- tempUMAP2$State
        cat("\n cells ordered ...\n")
        
        if( ( opt$no_of_paths != "none" ) & ( opt$start_cluster != "none" ) ){
            cat("\nRe-Ordering cells based on desired starting cluster and number of possible paths ...\n")
            #pData(tempUMAP2)$State <- factor(tempUMAP2[[opt$metadata_use]])
            tempUMAP2 <- orderCells(tempUMAP2 , root_state = opt$start_cluster,num_paths = as.numeric(opt$no_of_paths))
            
        } else if( opt$no_of_paths != "none" ){
            cat("\nRe-Ordering cells based on number of possible paths ...\n")
            tempUMAP2 <- orderCells(tempUMAP2 , num_paths = as.numeric(opt$no_of_paths))
            
        } else if( opt$start_cluster != "none" ){
            cat("\nRe-Ordering cells based on desired starting cluster ...\n")
            #tempUMAP2$State <- as.factor(tempUMAP2[[opt$metadata_use]])
            tempUMAP2 <- orderCells(tempUMAP2 , root_state = opt$start_cluster)
        }
        
        #tempUMAP2$State <- orig.state
        
        cat("\nReassigning gene expression to object \n")
        tempUMAP3 <- tempUMAP2
        tempUMAP3@assayData$exprs <- temp@assayData$exprs
        tempUMAP3@featureData <- temp@featureData
        tempUMAP3 <- estimateSizeFactors(tempUMAP3)
        
        saveRDS(tempUMAP3,paste0(opt$output_path,"/",reduc,"/Monocle_trajectory_object.rds"))
    }
    #---------
    
    
    
    ########################
    ### PLOTTING RESULTS ###
    ########################
    cat("\nPlotting results ...\n")
    
    for(i in c(unlist(strsplit(opt$metadata_use,",")),"State","Pseudotime") ){
        message(i)
        pdf(paste0(opt$output_path,"/",reduc,"/Trajectories_by_",i,".pdf"),width = 8,height = 8,useDingbats = F)
        print(plot_cell_trajectory(tempUMAP3, color_by = i,y=2,x=1,cell_size = .5))
        
        print(plot_complex_cell_trajectory(tempUMAP3, color_by = i, show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3))
        
        x <- factor(tempUMAP3@phenoData@data[,i])
        try(scatter3D(tempUMAP3@reducedDimS[1,],tempUMAP3@reducedDimS[2,],tempUMAP3@reducedDimS[3,], bty = "g", theta = 40,
        colvar = as.integer(factor(tempUMAP3@phenoData@data[,i])), col =  hue_pal()(length(levels(x))),
        pch=16,cex=.5,xlab = "Dim1", ylab ="Dim2", zlab = "Dim3",colkey = list(at = seq(1.5,length(levels(x))+.5,length.out = length(levels(x))),
        side = 1, labels = levels(x) ) ) )
        try(scatter3D(tempUMAP3@reducedDimK[1,],tempUMAP3@reducedDimK[2,],tempUMAP3@reducedDimK[3,],add=T,col="black",pch=16,cex=.5))
        
        for(j in names(DATA@reductions) ){
            message("   ",j)
            plot(DATA@reductions[[j]]@cell.embeddings[,1:2],col=hue_pal()(length(levels(x)))[x],cex=.7,pch=16,main=i,xlab = paste0(j,"_1"),ylab = paste0(j,"_2"),las=1)
        }
        
        dev.off()
    }
    #---------
    
    
    
    if( casefold(opt$diff_testing) %in% c("true","yes") ){
        
        
        #########################################
        ### COMPUTING DGE PER BRANCHING POINT ###
        #########################################
        cat("\nComputing DGE per Branch ...\n")
        n_branches <- length(tempUMAP3@auxOrderingData$DDRTree$branch_points)
        
        for(i in 1:n_branches){
            message('processing branch ',i)
            if(file.exists(paste0(opt$output_path,"/",reduc,"/diff_test_branch_",i,".csv"))){
                BEAM_res <- read.csv2(paste0(opt$output_path,"/",reduc,"/diff_test_branch_",i,".csv"),row.names = 1)
            } else {
                BEAM_res <- BEAM(tempUMAP3, branch_point = i, cores = detectCores()-1)
                BEAM_res <- BEAM_res[order(BEAM_res$qval),]
                write.csv2(BEAM_res,paste0(opt$output_path,"/",reduc,"/diff_test_branch_",i,".csv"),row.names = T)
            }
            n <- min(100, length(row.names(subset(BEAM_res, qval < 1e-6))))
            pdf(paste0(opt$output_path,"/",reduc,"/Trajectory_branch",i,".pdf"),width = 5,height = 7,useDingbats = F)
            try(plot_genes_branched_heatmap(tempUMAP3[rownames(BEAM_res)[1:n],],branch_point = i, norm_method= "log",
            num_clusters = 4, use_gene_short_name = T,  show_rownames = T))
            dev.off()
        }
        
        #---------
        
    }
    
    
    ###################################
    ### SAVING RAW Seurat.v3 OBJECT ###
    ###################################
    cat("\n### Saving the updated back to the Seurat object used as input###\n")
    DATA$Pseudotime <- tempUMAP3$Pseudotime
    DATA$State <- tempUMAP3$State
    
    #saveRDS(DATA, file = opt$Seurat_object_path )
    write.csv2(DATA@meta.data,paste0(opt$output_path,"/",reduc,"/Metadata_with_trajectory.csv"),row.names = T)
    #---------
    
}




if( opt$method_use == 'slingshot' ){
    
    library(Seurat)
    library(destiny)
    library(slingshot)
    library(scater)
    library(scales)
    library(mclust)
    library(RColorBrewer)
    a <- readRDS("/Users/paulo.barenco/Desktop/NBIS/Projects/M_Kasper_1709/trajectory_analysis_epithelial_190917/analysis/integrated_1000_k10/clustering/Seurat_object.rds")
    #a <- a[,!(a$louvain_0.9 %in% c(3,5)) ]
    dim(a)
    
    dm <- DiffusionMap(t(a@assays$mnn@data),k = 10,n_eigs = 15)
    a@reductions[["dm"]] <- CreateDimReducObject(embeddings = as.matrix(dm@eigenvectors),key = "DC_",assay = "RNA")
    
    k <- 17
    cl1 <- factor(Mclust(a@reductions[["umap10"]]@cell.embeddings,G = 1:k)$classification)
    a@meta.data[[paste0("mclust",k)]] <- cl1
    DimPlot(a,reduction = "umap",group.by = paste0("mclust",k),label = T)
    DimPlot(a,reduction = "umap10",group.by = "louvain_0.9",label = T)
    plot(a@reductions[["dm"]]@cell.embeddings[,1:2], col = scales::hue_pal()(length(unique(cl1)))[cl1], asp = 1, pch = 16,cex=.5)
    plot(a@reductions[["dm"]]@cell.embeddings[,3:4], col = scales::hue_pal()(length(unique(cl1)))[cl1], asp = 1, pch = 16,cex=.5)
    plot(a@reductions[["dm"]]@cell.embeddings[,5:6], col = scales::hue_pal()(length(unique(cl1)))[cl1], asp = 1, pch = 16,cex=.5)
    plot(a@reductions[["dm"]]@cell.embeddings[,7:8], col = scales::hue_pal()(length(unique(cl1)))[cl1], asp = 1, pch = 16,cex=.5)
    plot(a@reductions[["dm"]]@cell.embeddings[,9:10], col = scales::hue_pal()(length(unique(cl1)))[cl1], asp = 1, pch = 16,cex=.5)
    
    lin1 <- getLineages(a@reductions[["umap10"]]@cell.embeddings, cl1,start.clus=2)
    
    #plot curve centroids
    red <- a@reductions[["umap"]]@cell.embeddings
    centroids <-  sapply( unique(cl1) , red=red, cl1=cl1, function(x,red,cl1) { colMedians( red[cl1==x,] ) })
    colnames(centroids) <- as.character(unique(cl1))
    
    plot(a@reductions[["umap"]]@cell.embeddings, col = scales::hue_pal()(length(unique(cl1)))[cl1], asp = 1, pch = 16,cex=.5)
    for(i in lin1@lineages){  lines( t(centroids)[ i, ] , lwd=1.5) }
    points(t(centroids)[,1:2],bg="#ffffff95",cex=2,pch=21)
    text(t(centroids)[,1:2],labels = colnames(centroids),cex=.6)
    
    
    
    
    

        
    
    
    
    
    
    
    
    
    colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
    plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
    plot(a@reductions[["umap"]]@cell.embeddings, col = scales::hue_pal()(length(unique(cl1)))[cl1], asp = 1, pch = 16,cex=.5)
    lines(crv1, lwd = 3, col = 'black')
    
    
    
    
    
    
    sce <- slingshot(sim, clusterLabels = "louvain_0.9", reducedDim = 'UMAP10',
    start.clus= as.character(7), end.clus = as.character(c(6,8,4,0)))
    
    sce
    plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
    lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')
    
    
    
    rd <- a@reductions[["umap10"]]@cell.embeddings
    k <- 12
    cl1 <- factor(Mclust(rd,G = 1:k)$classification)
    #cl1 <- factor(kmeans(rd, centers = 20)$cluster)
    #cl1 <- factor(a$louvain_0.9)
    a@meta.data[[paste0("mclust",k)]] <- cl1
    # plot(a@reductions[["umap"]]@cell.embeddings, col = scales::hue_pal()(length(unique(cl1)))[cl1], pch=16, asp = 1)
    
    plot(rd, col = scales::hue_pal()(length(unique(cl1)))[cl1], pch=16, asp = 1)
    DimPlot(a,reduction = "umap",group.by = paste0("mclust",k),label = T)
    #cl <- factor(a$res.1.15)
    
    
    rd <- a@reductions[["umap10"]]@cell.embeddings
    cl1 <- a@meta.data[[paste0("mclust",k)]]
    
    lin1 <- getLineages(a@reductions[["umap10"]]@cell.embeddings, cl1,
    start.clus= as.character(2), end.clus = as.character(c(11,9,4,8,10,12)))
    lin1@reducedDim <- a@reductions[["umap"]]@cell.embeddings
    crv1 <- getCurves(lin1,stretch=1)
    
    
    
    # plot(rd, col = scales::hue_pal()(length(unique(cl1)))[cl1], asp = 1, pch = 16,cex=.5)
    # lines(lin1, lwd = 1, col = 'black')
    
    plot(a@reductions[["umap"]]@cell.embeddings, col = scales::hue_pal()(length(unique(cl1)))[cl1], asp = 1, pch = 16,cex=.5)
    lines(crv1, lwd = 3, col = 'black')
    
    
    plot(a@reductions[["umap"]]@cell.embeddings, col = colorRampPalette(c("grey","navy"))(ncol(a)+2)[crv1@curves$curve1$ord], asp = 1, pch = 16,cex=.5)
    lines(crv1, lwd = 3, col = 'black')
    
    
    
    pheatmap::pheatmap(log2(lin1@slingParams$dist+1) )
    
    plot(rd, col = scales::hue_pal()(length(unique(cl1)))[cl1], asp = 1, pch = 16,cex=.5)
    lines(lin1, lwd = 1, col = 'black')
    
    
    
    
    #Test branching point 1 (clusters 1-2 vs 7-2)
    root <- 1
    leafs <- c(2,7)
    
    crv1@lineages
    to_test <- cl1 %in% c(1,2,7)
    Y <- a@assays$RNA@data[,to_test]
    Y <- Y[rowSums(Y>0)>100,]
    dim(Y)
    Group <- factor( paste0( "clust",cl1[cl1 %in% c(1,2,7)]) )
    design <- model.matrix( ~ Group)
    colnames(design)
    contr <- c( -1,1,1)
    
    #Compute a "Blocking" experimental design to find any differences between strains conserved among tissues
    y <- DGEList(counts=Y, group=Group)
    
    keep <- rowSums(cpm(y)>1) >= 3        #Here, we need to filter out genes containing 0 reads again (since we filtered the samples)
    y <- y[keep, , keep.lib.sizes=FALSE]
    
    y <- calcNormFactors(y, method = "TMM")
    y <- estimateDisp(y, design, robust=TRUE)
    
    fit <- glmQLFit(y, design, robust = T) #Better for low number of samples or 1 sample
    normcounts <- cpm(y, normalized.lib.sizes = T, log = T)
    lrt1 <- glmQLFTest(fit,coef = 2)
    DEtable <- as.matrix(topTags(lrt1,n = "all")$table)
    
    lrt2 <- glmQLFTest(fit,coef = 3)
    DEtable2 <- as.matrix(topTags(lrt2,n = "all")$table)
    
    
    top <- DEtable[1:200,"logFC"] - DEtable2[1:200,"logFC"]
    top <- names(sort(abs(top)))[1:100]
    
    pheatmap( cbind( Y[rownames(DEtable)[1:100],order(Group)] ),scale = "row" )
    
    
    
    
    cat("\nComputing DGE per Branch ...\n")
    n_branches <- length(tempUMAP3@auxOrderingData$DDRTree$branch_points)
    
    for(i in 1:n_branches){
        message('processing branch ',i)
        if(file.exists(paste0(opt$output_path,"/",reduc,"/diff_test_branch_",i,".csv"))){
            BEAM_res <- read.csv2(paste0(opt$output_path,"/",reduc,"/diff_test_branch_",i,".csv"),row.names = 1)
        } else {
            BEAM_res <- BEAM(tempUMAP3, branch_point = i, cores = detectCores()-1)
            BEAM_res <- BEAM_res[order(BEAM_res$qval),]
            write.csv2(BEAM_res,paste0(opt$output_path,"/",reduc,"/diff_test_branch_",i,".csv"),row.names = T)
        }
        n <- min(100, length(row.names(subset(BEAM_res, qval < 1e-6))))
        pdf(paste0(opt$output_path,"/",reduc,"/Trajectory_branch",i,".pdf"),width = 5,height = 6,useDingbats = F)
        try(plot_genes_branched_heatmap(tempUMAP3[rownames(BEAM_res)[1:n],],branch_point = i, norm_method= "log",
        num_clusters = 4, use_gene_short_name = T,  show_rownames = T))
        dev.off()
    }
    
    
    
    
    crv1@slingParams$
    cl1 %in% crv1@lineages$Lineage1
    crv1@lineages$Lineage1
    
    crv1@curves$curve1$ord
    crv1@curves$curve2$ord
    
    t <- log2(crv1@curves$curve1$ord)
    t[!cl1 %in% crv1@lineages$Lineage1] <- NA
    
    # for time, only look at the 100 most variable genes
    Y <- a@assays$RNA@data[a@assays$RNA@var.features,]
    
    # fit a GAM with a loess term for pseudotime
    require(gam)
    require(clusterExperiment)
    gam.pval <- apply(Y,1,t=t,function(z,t){
        d <- data.frame(z=z, t=t)
        tmp <- gam(z ~ lo(t), data=d)
        p <- summary(tmp)[3][[1]][2,3]
        p
    })
    topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
    heatdata <- as.matrix(Y[topgenes, order(t, na.last = NA)])
    heatclus <- cl1[order(t, na.last = NA)]
    ce <- ClusterExperiment(heatdata, heatclus)
    png(paste0(opt$output_path,"/heatmap_DEG2.png"),width = 1000,height = 1500,res=200 )
    plotHeatmap(ce, clusterSamplesData = "orderSamplesValue",
    visualizeData = 'original')
    dev.off()
    
    
    
    
    
    
    
    
    
    
    
    
}



#############################
### SYSTEM & SESSION INFO ###
#############################
#---------
print_session_info()
#---------
