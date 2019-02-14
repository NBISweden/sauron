#!/usr/bin/env Rscript

### LOAD OPTPARSE
#---------
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')};
library(optparse)
#---------


### DEFINE PATH TO LOCAL FILES
#---------
cat("\nRunnign DIMENS> REDUCTION AND CLUSTERING with the following parameters ...\n")
option_list = list(
  make_option(c("-i", "--objects_paths"),         type = "character",   metavar="character",   default='none',  help="Path to the Seurat object"),
  make_option(c("-n", "--object_names"),          type = "character",   metavar="character",   default='none',  help="Column names in the Metadata matrix (only factors allowed, not continuous variables)"),
  make_option(c("-c", "--object_clusters"),       type = "character",   metavar="character",   default='none',  help="Variables to be regressed out using linear modeling."),
  make_option(c("-m", "--metadata_use"),          type = "character",   metavar="character",   default='none',  help="Variables to be regressed out using linear modeling."),
  make_option(c("-s", "--species_use"),          type = "character",   metavar="character",   default='none',  help="Variables to be regressed out using linear modeling."),
  make_option(c("-d", "--lig_recp_database"),     type = "character",   metavar="character",   default='none',  help="Batch-correction method to be used. 'MNN', 'Scale' and 'Combat' are available at the moment. The batches (column names in the metadata matrix) to be removed should be provided as arguments comma separated. E.g.: 'Combat,sampling_day'. For MNN, an additional integer parameter is supplied as the k-nearest neighbour."),
  make_option(c("-l", "--ligand_objects"),        type = "character",   metavar="character",   default='none',  help="Batch-correction method to be used. 'MNN', 'Scale' and 'Combat' are available at the moment. The batches (column names in the metadata matrix) to be removed should be provided as arguments comma separated. E.g.: 'Combat,sampling_day'. For MNN, an additional integer parameter is supplied as the k-nearest neighbour."),
  make_option(c("-r", "--receptor_objects"),      type = "character",   metavar="character",   default='top,5', help="Method and threshold level for selection of significant principal components. The method should be separated from the threshold via a comma. 'top,5' will use the top 5 PCs, which is the default. 'var,1' will use all PCs with variance above 1%."),
  make_option(c("-f", "--aux_functions_path"),    type = "character",   metavar="character",   default='none',  help="File with supplementary functions"),
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
#---------





### LOAD LIBRARIES
#---------
cat("\nLoading/installing libraries ...\n")
source(opt$aux_functions_path)
pkgs <- c("Seurat","rafalib","scran","biomaRt","scater","dplyr","RColorBrewer","dbscan","flowPeaks","scales","igraph")
inst_packages(pkgs)
#---------





### LOAD Seurat OBJECTS
#---------
cat("\nLoading/ data ...\n")
objects_paths <- unlist(strsplit(opt$objects_paths,","))
object_names <- unlist(strsplit(opt$object_names,","))

print(objects_paths)
print(object_names)

for(i in 1:length(object_names)){
  assign( object_names[i] , readRDS(objects_paths[i]) )
}
#---------





### Filtering clusters for each dataset
#---------
cat("\nThe following clusters will be used for the respective datasets ...\n")
object_clusters <- unlist(strsplit(opt$object_clusters,";"))
object_clusters <- lapply( object_clusters , function(x) unlist(strsplit(x,",")) )
names(object_clusters) <- object_names
print(object_clusters)

for(i in object_names){
  if( length(object_clusters[[i]]) >= 1){
    #Assign cell identities based on the metadata
    temp <- get(i)
    temp@ident <- factor(NULL)
    temp <- SetIdent(temp,ident.use = temp@meta.data[ , object_clusters[[i]][1] ])

    #Filter cells from a clusters
    if( length(object_clusters[[i]]) >= 2){
      head(temp@meta.data[ , object_clusters[[i]][1] ])
      head(object_clusters[[i]][-1])
      cell_use <- rownames(temp@meta.data) [ temp@meta.data[ , object_clusters[[i]][1] ] %in% object_clusters[[i]][-1] ]
      temp <- SubsetData(temp,cells.use = cell_use)
      dim(temp@meta.data)
      dim(temp@data)
    }
    
    assign( object_names[i] , temp )
    
  }
}
rm(temp)
gc()
#---------


#Import receptor-ligand interaction dataset
#---------
cat("\nLoading/ receptor-ligand interaction dataset ...\n")

#import table
L_R_pairs <- read.csv2("/Users/paulo.barenco/Box/repos/single_cell_analysis/support_files/ligand_receptor/ligand_receptor_pairs.csv")

#convert symbols to mouse mgi IDs
if(opt$species_use == "mouse"){ human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")    ;    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
m_h_symbols = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = unique(c(as.character(L_R_pairs$ligand),as.character(L_R_pairs$receptor))) , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=F,valuesL = "hgnc_symbol") 

#replace them in the matrix and remove NAs
L_R_pairs$ligand = m_h_symbols[match(L_R_pairs$ligand, m_h_symbols[,2]),1]
L_R_pairs$receptor = m_h_symbols[match(L_R_pairs$receptor, m_h_symbols[,2]),1]
}


L_R_pairs <- na.omit(L_R_pairs)
print(head(L_R_pairs))
dim(L_R_pairs)

#png(filename = paste0(output_path,"/VENN_Fibroblasts_Keratinocyte_ligand_receptors.png"),width = 1300,height = 1000,res = 150)
#mypar(1,2)

#VENN diagram for list of LIGANDS/RECEPTORS detected in Fibroblasts
#input  <- list("det_Fb"=det_Fb,"Lig"=L_R_pairs$Ligand.ApprovedSymbol,"Rec"=L_R_pairs$Receptor.ApprovedSymbol)
#a <- attr(venn::venn(input,zcolor = 1:3,lwd=3,col=1:3,cexil = c(1,1,1,1,1,1,2)*.9,cexsn = 1), "intersections")
#title(main = paste0("Venn diagram for\nLigand-Receptors detected\n in Fibroblasts"),line = -5)
#write.csv2(a$`det_Fb:Lig`,"Fibroblast_ligands.csv")
#write.csv2(a$`det_Fb:Rec`,"Fibroblast_receptors.csv")

#VENN diagram for list of LIGANDS/RECEPTORS detected in Fibroblasts
#input  <- list("det_Ep"=det_Ep,"Lig"=L_R_pairs$Ligand.ApprovedSymbol,"Rec"=L_R_pairs$Receptor.ApprovedSymbol)
#b <- attr(venn::venn(input,zcolor = 1:4,lwd=3,col=1:4,cexil = c(1,1,1,1,1,2,1)*.9,cexsn = 1), "intersections")
#title(main = paste0("Venn diagram for\nLigand-Receptors detected\n in Keratinocytes"),line = -5)
#write.csv2(b$`det_Ep:Lig`,"Keratinocyte_ligands.csv")
#write.csv2(b$`det_Ep:Rec`,"Keratinocyte_receptors.csv")
#dev.off()

#det_Fb_lig <- L_R_pairs$Ligand.ApprovedSymbol %in% a$`det_Fb:Lig`
#det_Ep_rec <- L_R_pairs$Receptor.ApprovedSymbol %in% b$`det_Ep:Rec`
#Fib_L_R_pairs <- L_R_pairs[det_Fb_lig & det_Ep_rec,]
#write.csv2(Fib_L_R_pairs,"Fibroblast(L)-keratinocyte(R)_pairs.csv")

#det_Ep_lig <- L_R_pairs$Ligand.ApprovedSymbol %in% b$`det_Ep:Lig`
#det_Fb_rec <- L_R_pairs$Receptor.ApprovedSymbol %in% a$`det_Fb:Rec`
#Ep_L_R_pairs <- L_R_pairs[det_Ep_lig & det_Fb_rec,]
#write.csv2(Ep_L_R_pairs,"Keratinocyte(L)-fibroblast(R)_pairs.csv")
#---------







### Identifying relevant markers across embrionic age for each EPITHELIAL population
#---------
cat("\nIdentifying detected genes in each cluster ...\n")

for(m in object_names){
  message(paste0("\n\nProcessing dataset ",m," ..."))
  
  #Create empty objects to store results
  marker_list <- data.frame()
  cluster_data <- list()
  
  #if(!dir.exists(paste0(output_path,"/DGE_embrionic.age_",m))){dir.create(paste0(output_path,"/DGE_embrionic.age_",m))}
  
  for(i in unique(get(m)@ident)){    cat(paste0("... Processing cell cluster #",i," ..."),"\n")
    
    #Load only the cells for a specific cluster "i".
    temp <- SubsetData(get(m), cells.use = get(m)@cell.names[get(m)@ident == i]) #Select cell from a cluster
    temp@ident <- factor(NULL)
    #temp <- SetIdent(temp,ident.use = temp@meta.data$Embryonic_age)
    #temp@meta.data$Embryonic_age
    
    
    #Filter out non-detected genes
    det <- names(which(rowSums(as.matrix(temp@data) > 0.1) > ncol(temp@data)*.3  ))
    #if(m == "EPI"){use <- det[det %in% Fib_L_R_pairs$Receptor.ApprovedSymbol]
    #} else { use <- det[det %in% Fib_L_R_pairs$Ligand.ApprovedSymbol] }
    use <- det[det %in% c(L_R_pairs$ligand , L_R_pairs$receptor) ]
    temp_markers <- data.frame(gene=use,cluster=i)
    k <- as.matrix(temp@data)[use,]
    
    
    #Correlate expression with the time points as well as differences in expresssion due to cell numbers
    #temp_markers$cor.r <- sapply(as.character(temp_markers$gene),function(x)cor( c(k[x,],0,0,0), c(as.numeric(temp@meta.data$Embryonic_age),1,2,3), method = "pearson"))
    #temp_markers$cor.pvalue <- sapply(as.character(temp_markers$gene),function(x)cor.test(c(k[x,],0,0,0), c(as.numeric(temp@meta.data$Embryonic_age),1,2,3),method="pearson")$p.value)
    #temp_markers[is.na(temp_markers)] <- 0
    
    
    #Filter genes that have no trend in differential expression and neither belong to a cluster without significant change in cell number
    #temp_markers$signif_exp <- (temp_markers$cor.pvalue < 0.05) & (abs(temp_markers$cor.r) > 0.3)
    
    
    #marker_list[[i]] <- temp_markers
    marker_list <- rbind(marker_list, temp_markers)
    cluster_data[[i]] <- k[temp_markers$gene,]
    
    write.csv2(temp_markers,paste0(output_path,"/DGE_embrionic.age_",m,"/markers_embryonic.age_for_cluster",i,".csv"),row.names = T)
    png(filename = paste0(output_path,"/DGE_embrionic.age_",m,"/vioplot_for_cluster",i,".png"),width = 300*5,height = 300*1.5*length(unique(temp_markers$gene))/4,res = 150)
    print(VlnPlot(object = temp, features.plot = unique(temp_markers$gene), point.size.use = .1,nCol = 4)) #it does not work if you don't have the print command in front of it!
    dev.off()
  }
  
  marker_list$edge <- paste0(m,marker_list$cluster,">",marker_list$gene)
  #filt_marker_list <- marker_list[  rowSums(marker_list[,c("signif_exp","signif_cell","signif_sum"),]) > 0 ,]
  
  assign(paste0("marker_list_",m)  , marker_list)
  #assign(paste0("cors_",m)         , filt_marker_list)
  #assign(paste0("cluster_data_",m) , cluster_data)
  #rm(marker_list, filt_marker_list)
}
#---------



write.csv(get(paste0("marker_list_",m)) , paste0(opt$output_path,"/Marker_list_",m,".csv"),row.names = T)



