#!/usr/bin/env Rscript

### LOAD OPTPARSE
#---------
if(!require("optparse")){install.packages("optparse", repos='http://cran.us.r-project.org')};
library(optparse)
#---------


### DEFINE PATH TO LOCAL FILES
#---------
cat("\nRunnign LIGAND-RECEPTOR interactome script ...\n")
option_list = list(
  make_option(c("-i", "--objects_paths"),         type = "character",   metavar="character",   default='none',  help="Path to the Seurat object"),
  make_option(c("-n", "--object_names"),          type = "character",   metavar="character",   default='none',  help="Column names in the Metadata matrix (only factors allowed, not continuous variables)"),
  make_option(c("-c", "--object_clusters"),       type = "character",   metavar="character",   default='none',  help="List of clusters to be used."),
  make_option(c("-ml", "--metadata_ligands"),      type = "character",   metavar="character",   default='none',  help="Metadata to be used as parameter for scoring ligand edges."),
  make_option(c("-mr", "--metadata_receptors"),     type = "character",   metavar="character",   default='none',  help="Metadata to be used as parameter for scoring receptor edges."),
  make_option(c("-s", "--species_use"),           type = "character",   metavar="character",   default='none',  help="Specied to be used."),
  make_option(c("-d", "--lig_recp_database"),     type = "character",   metavar="character",   default='none',  help="Ligand-receptor matrix to look for gene pairs."),
  make_option(c("-l", "--ligand_objects"),        type = "character",   metavar="character",   default='none',  help="List of objects that will be used as ligands."),
  make_option(c("-r", "--receptor_objects"),      type = "character",   metavar="character",   default='top,5', help="List of objects that will be used as receptors."),
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
initial.options <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(sub("--file=","",initial.options[grep("--file=",initial.options)]))
source( paste0(script_path,"/inst_packages.R") )
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
  temp <- readRDS(objects_paths[i])
  temp@raw.data <- NULL
  assign( object_names[i] , temp )
}
rm(temp)
invisible(gc())
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
    
    
    #Filter cells from a clusters
    #if( length(object_clusters[[i]]) >= 2){
    #  head(temp@meta.data[ , object_clusters[[i]][1] ])
    #  head(object_clusters[[i]][-1])
    #  cell_use <- rownames(temp@meta.data) [ temp@meta.data[ , object_clusters[[i]][1] ] %in% object_clusters[[i]][-1] ]
    #  temp <- SubsetData(temp,cells.use = cell_use)
    #}
    
    temp@ident <- factor(NULL)
    temp <- SetIdent(temp,ident.use = temp@meta.data[ , object_clusters[[i]][1] ])
    assign( i , temp )
  }
}
rm(temp)
invisible(gc())
#---------


#Import receptor-ligand interaction dataset
#---------
cat("\nLoading/ receptor-ligand interaction dataset ...\n")

#import table
print(opt$lig_recp_database)
L_R_pairs <- read.csv2(opt$lig_recp_database,stringsAsFactors = F)

#convert symbols to mouse mgi IDs
if(opt$species_use == "mouse"){ human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")    ;    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
m_h_symbols = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = unique(c(as.character(L_R_pairs$ligand),as.character(L_R_pairs$receptor))) , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=F,valuesL = "hgnc_symbol") 

#replace them in the matrix and remove NAs
L_R_pairs$ligand = m_h_symbols[match(L_R_pairs$ligand, m_h_symbols[,2]),1]
L_R_pairs$receptor = m_h_symbols[match(L_R_pairs$receptor, m_h_symbols[,2]),1]
}

L_R_pairs$ligand <- as.character(L_R_pairs$ligand)
L_R_pairs$receptor <- as.character(L_R_pairs$receptor)


L_R_pairs <- na.omit(L_R_pairs)
write.csv(L_R_pairs , paste0(opt$output_path,"/Ligand_Receptor_pairs.csv"),row.names = T)


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
obj_list <- list(L = unlist(strsplit(opt$ligand_objects,",")),
                 R = unlist(strsplit(opt$receptor_objects,",")))

for( t in names(obj_list) ){
  message(paste0("\n\nProcessing ",t," ..."))
  #Create empty objects to store results
  marker_list <- data.frame()
  cluster_data <- list()
  
  for(m in obj_list[[t]]){
    message(paste0("\n\nProcessing dataset ",m," ..."))

    #if(!dir.exists(paste0(output_path,"/DGE_embrionic.age_",m))){dir.create(paste0(output_path,"/DGE_embrionic.age_",m))}
    temp <- get(m)
    temp@ident <- factor(NULL)
    temp <- SetIdent(temp,ident.use = temp@meta.data[ , object_clusters[[m]][1] ])
    if(length(object_clusters[[m]]) > 1){
      clust_use <- object_clusters[[m]][-1]
    } else { clust_use <- unique(temp@ident) }

    for(i in clust_use){    cat(paste0("... Processing cell cluster #",i," ..."),"\n")
      
      #Load only the cells for a specific cluster "i".
      temp2 <- SubsetData(temp, cells.use = temp@cell.names[temp@ident == i]) #Select cell from a cluster
      #temp <- SetIdent(temp,ident.use = temp@meta.data$Embryonic_age)
      #temp@meta.data$Embryonic_age

      #Filter out non-detected genes
      det <- rownames(temp2@data) [ rowSums(as.matrix(temp2@data) > 0) > ncol(temp2@data)*.5 ]
      #if(m == "EPI"){use <- det[det %in% Fib_L_R_pairs$Receptor.ApprovedSymbol]
      #} else { use <- det[det %in% Fib_L_R_pairs$Ligand.ApprovedSymbol] }
      
      if(  t == "L" ) { use <- det[det %in% L_R_pairs$ligand ]
      } else {  use <- det[det %in% L_R_pairs$receptor ] }
      
      temp_markers <- data.frame(gene=use,cluster=i)
      temp2@data <- temp2@data[use,]
  
      #Compute differential expression with the time points as well as differences in expresssion due to cell numbers
      
      if( sum(c(opt$metadata_ligands,opt$metadata_receptors) %in% colnames(temp2@meta.data)  ) > 0 ){
        k <- as.matrix(temp2@data)[use,]
        if(  t == "L" ) {   meta <- c(as.numeric(temp2@meta.data[,opt$metadata_ligands]))
        } else {  meta <- c(as.numeric(temp2@meta.data[,opt$metadata_receptors]))}

        cors <- lapply(as.character(temp_markers$gene),function(x)cor.test( c(k[x,]), meta, method="pearson"))
        temp_markers$cor.r <- sapply(cors,function(x) x$estimate )
        temp_markers$cor.pvalue <- sapply(cors,function(x) x$p.value )
        temp_markers[is.na(temp_markers)] <- 0
        temp_markers$signif_exp <- (temp_markers$cor.pvalue < 0.05) & (abs(temp_markers$cor.r) > 0.3)
      }
      
      if( t == "L" ) { temp_markers$edge <- paste0(t,"_",m,"_",temp_markers$cluster,">",temp_markers$gene)
      } else { temp_markers$edge <- paste0(temp_markers$gene,">",t,"_",m,"_",temp_markers$cluster) }
      
      #marker_list[[i]] <- temp_markers
      marker_list <- rbind(marker_list, temp_markers)
      print(dim(marker_list))
      #cluster_data[[i]] <- k[temp_markers$gene,]
      
      #write.csv2(temp_markers,paste0(opt$output_path,"/DGE_embrionic.age_",m,"/markers_embryonic.age_for_cluster",i,".csv"),row.names = T)
      #png(filename = paste0(opt$output_path,"/vioplot_for_cluster",i,".png"),width = 300*5,height = 300*1.5*length(unique(temp_markers$gene))/4,res = 150)
      #print(VlnPlot(object = temp2, features.plot = unique(temp_markers$gene), point.size.use = .1,nCol = 4)) #it does not work if you don't have the print command in front of it!
      #dev.off()
    }
  }
  
  cat("\nAssigning and saving results to list ...\n")
  #filt_marker_list <- marker_list[  rowSums(marker_list[,c("signif_exp"),]) > 0 ,]
  
  #compute uniqueness
  marker_list$uniqueness <- sapply( marker_list$gene, function(x) sum(marker_list$gene == x) )
  marker_list$uniqueness <- marker_list$uniqueness
  
  assign(paste0("marker_list_",t)  , marker_list)
  #assign(paste0("cors_",m)         , filt_marker_list)
  #assign(paste0("cluster_data_",m) , cluster_data)
  #rm(marker_list, filt_marker_list)
  
  write.csv2(get(paste0("marker_list_",t)) , paste0(opt$output_path,"/Marker_list_",t,".csv"),row.names = T)
}
marker_list_all <- rbind(marker_list_L,marker_list_R)
rm(temp,temp2)
invisible(gc())
#---------

#marker_list_L <- read.csv2('~/Downloads/filtered_gene_bc_matrices/analysis/5-interactome/Marker_list_L.csv',row.names = 1,stringsAsFactors = F)
#marker_list_R <- read.csv2('~/Downloads/filtered_gene_bc_matrices/analysis/5-interactome/Marker_list_R.csv',row.names = 1,stringsAsFactors = F)



#CALCULATING THE CONNECTIONS BETWEEEN CELLS THAT SHARE LIGAND-RECEPTOR INTERACTIONS
#---------
cat("\nCalculating connection between ligands and receptors ...\n")
k <- data.frame(t(as.data.frame(strsplit( c(marker_list_L$edge, marker_list_R$edge),split = ">"))))
all_genes <- unique(c(L_R_pairs$ligand,L_R_pairs$receptor))
datasets <- as.character(unique(unlist(k)[! (unlist(k) %in% all_genes)]))

colnames(k) <- c('ligand','receptor')
#k <- k[c(cors_FIB$signif_exp, cors_EPI$signif_cell),]
#k <- data.frame(L=c(L_R_pairs$ligand,k[,1]) , R=c(L_R_pairs$receptor,k[,2]),stringsAsFactors = F)
det_L_R_pairs <- L_R_pairs[ (L_R_pairs$ligand %in% as.character(k[,2])) & (L_R_pairs$receptor %in% as.character(k[,1])) , ]
k <- rbind(det_L_R_pairs , k)
rownames(det_L_R_pairs) <- 1:nrow(det_L_R_pairs)
#---------




#GRAPH_1 - All
#---------
cat("\nPlotting raw graph ...\n")
#Plot the Graph network of interacting Ligand-Receptors for fibroblasts and Epithelial cells
all_edges <- unlist(strsplit(paste(k$ligand,k$receptor,collapse = " ")," "))
g <- graph( edges=all_edges, directed=T )
g <- simplify(g, remove.multiple = T, remove.loops = T)

cols <- ifelse( V(g)$name %in% L_R_pairs$ligand,hue_pal()(10)[8],
                ifelse( V(g)$name %in% L_R_pairs$receptor,hue_pal()(10)[4], hue_pal()(10)[1]) )

png(filename = paste0(opt$output_path,"/Graph_all_ligand_receptor.png"),width = 1200,height =1200,res = 100)
plot.igraph(g, vertex.label.color="black",vertex.label.family="sans",vertex.label.cex=1,vertex.label.font=2,
            vertex.shape="circle", vertex.size=10,edge.arrow.width=1,edge.arrow.size=.5,
            vertex.color=paste0(cols,90 ),
            vertex.frame.color=cols )
dev.off()
#---------




#GRAPH_2 - Filtered & colored by uniqueness
#---------
cat("\nPlotting filtered graph ...\n")
#Plot the Graph network of interacting Ligand-Receptors for fibroblasts and Epithelial cells
k <- k[ (k$ligand %in% k$receptor) | (k$receptor %in% k$ligand) , ]
all_edges <- unlist(strsplit(paste(k$ligand,k$receptor,collapse = " ")," "))
g <- graph( edges=all_edges, directed=T )
g <- simplify(g, remove.multiple = T, remove.loops = T)
cols <- ifelse( V(g)$name %in% L_R_pairs$ligand,hue_pal()(10)[8],
                ifelse( V(g)$name %in% L_R_pairs$receptor,hue_pal()(10)[4], hue_pal()(10)[1]) )
l <- layout_with_sugiyama(g)
l <- l$layout[,2:1]
for (i in unique(l[,1])){  r <- rank( l[l[,1] == i ,2 ] ) ; l[l[,1] == i ,2 ] <- r / ( max(r)+1 ) }

print(length(unique(unlist(k))))

edge_info <- rbind(marker_list_L,marker_list_R)[match(paste(as_edgelist(g)[,1],as_edgelist(g)[,2],sep = ">"), rbind(marker_list_L,marker_list_R)$edge),"uniqueness"]
edge_color <- ifelse(is.na(edge_info),"black",ifelse(edge_info == 1,"black",ifelse(edge_info == 2,"grey50","grey70")))
edge_scalar <- ifelse(is.na(edge_info),.5,ifelse(edge_info == 1,2,ifelse(edge_info == 2,2,.5)))

png(filename = paste0(opt$output_path,"/Graph_paired_ligand_receptor_uniqueness.png"),width = 1800,height = 1800*max(.6, length(unique(unlist(k)))/100 ),res = 150)
plot.igraph(g, vertex.label.color="black",vertex.label.family="sans",vertex.label.cex=1,vertex.label.font=2,
            vertex.shape="vrectangle", vertex.size=25,vertex.size2=4/max(.6, length(unique(unlist(k)))/100 ),edge.arrow.width=edge_scalar*4,edge.arrow.size=0,
            vertex.color=paste0(cols,90 ), vertex.frame.color=cols ,layout=-l,
            edge.color=edge_color,edge.width=edge_scalar,asp = max(.6, length(unique(unlist(k)))/100 )  )
legend(-1,1.4,legend = c("clusters","ligands","receptors"),pch=22,xjust = 0,yjust = 1,
       pt.bg = paste0(hue_pal()(10)[c(1,8,4)],90),bty = "n", col=hue_pal()(10)[c(1,8,4)] )
legend(-.5,1.4,legend = c("specific* (1 connection)","medium* (2 connections)","broad* (>2 connections)"),lty=1, xjust = 0,yjust = 1,lwd=c(2,2,1),
       col = c("black","grey50","grey70"),bty = "n")
text(-1,-1.2,"* only the edges between genes and clusters were weighted.",adj = 0)
dev.off()
#---------







#GRAPH_3 - Filtered & colored by correlation to metadata
#---------
# cor_pal <- colorRampPalette(c("blue","navy","grey80","firebrick3","red"))(19)
# myedges <- apply(as_edgelist(g),1,function(x) paste(x,collapse = ">"))
# edge_cors <- marker_list_all[match(myedges,marker_list_all$edge),"cor.r"]
# mycolor <- ifelse( is.na(edge_cors) , "black" , cor_pal[round( (edge_cors+1)*9+1,0)] )
# 
# png(filename = paste0(opt$output_path,"/Graph_paired_ligand_receptor_metadata.png"),width = 1800,height = 1800*max(.6, length(unique(unlist(k)))/100 ),res = 150)
# plot.igraph(g, vertex.label.color="black",vertex.label.family="sans",vertex.label.cex=1,vertex.label.font=2,
#             vertex.shape="vrectangle", vertex.size=25,vertex.size2=4/max(.6, length(unique(unlist(k)))/100 ),edge.arrow.width=edge_scalar*4,edge.arrow.size=0,
#             vertex.color=paste0(cols,90 ), vertex.frame.color=cols ,layout=-l,
#             edge.color=mycolor,edge.width=3,asp = max(.6, length(unique(unlist(k)))/100 )  )
# legend(-1,1.4,legend = c("clusters","ligands","receptors"),pch=22,xjust = 0,yjust = 1,
#        pt.bg = paste0(hue_pal()(10)[c(1,8,4)],90),bty = "n", col=hue_pal()(10)[c(1,8,4)] )
# legend(-.5,1.4,legend = c("specific* (1 connection)","medium* (2 connections)","broad* (>2 connections)"),lty=1, xjust = 0,yjust = 1,lwd=c(2,2,1),
#        col = c("black","grey50","grey70"),bty = "n")
# text(-1,-1.2,"* only the edges between genes and clusters were weighted.",adj = 0)
# dev.off()
#---------






#Calculating all possible paths between clusters
#---------
cat("\nCalculating all possible paths between clusters ...\n")
ligand_clusters <- V(g)$name[grepl("L_",V(g)$name)]
receptor_clusters <- V(g)$name[grepl("R_",V(g)$name)]
all_paths <- data.frame()

for(i in ligand_clusters){
  temp <- all_simple_paths(g, from = i, to = receptor_clusters)
  temp <- lapply(temp,function(x) x$name)
  temp <- t(as.data.frame(temp))
  all_paths <- rbind(all_paths, temp)
}

colnames(all_paths) <- c( "ligand_clusters", "ligand" , "receptor" , "receptor_clusters")
rownames(all_paths) <- 1:nrow(all_paths)
write.csv2(all_paths , paste0(opt$output_path,"/All_paths.csv"),row.names = T)
#---------






#Summarizing all possible paths between clusters
#---------
cat("\nSummarizing all possible paths between clusters ...\n")

path_summary <- lapply(unique(all_paths[,'ligand_clusters']), function(x) c(table(all_paths[all_paths[,'ligand_clusters']==x,'receptor_clusters'])) )
names(path_summary) <- unique(all_paths[,'ligand_clusters'])  ;   path_summary <- t(as.matrix(as.data.frame(path_summary)))
print(path_summary)

res <- data.frame()
for(i in unique(all_paths[,'receptor_clusters']) ){
  for(j in unique(all_paths[,'ligand_clusters']) ){
    temp <- t(data.frame(setNames(c( i , j , path_summary[j,i] ),c("rec","lig","inter"))))
    res <- rbind(res , temp,deparse.level = 0) }}
rownames(res) <- 1:nrow(res)


g2 <- graph_from_data_frame(res)
l <- layout_with_sugiyama(g2)  ;  l <- scale(-l$layout[,2:1])
cols <- ifelse( V(g2)$name %in% unique(all_paths[,'ligand_clusters']),hue_pal()(10)[8],hue_pal()(10)[4])

png(filename = paste0(opt$output_path,"/Interaction_count.png"),width = 1200,height = 800,res = 150)
plot.igraph(g2, vertex.label.color="black",vertex.label.family="sans",vertex.label.cex=1,vertex.label.font=2,
            vertex.shape="vrectangle", vertex.size=50,vertex.size2=15,edge.arrow.size=0,
            vertex.color=paste0(cols,90 ), vertex.frame.color=cols ,layout=-l,
            edge.width=as.numeric(as.character(res[,3]))/max(as.numeric(as.character(res[,3])))*6,
            edge.label=res[,3],edge.label.font=2,edge.label.cex=1,asp=.6)
dev.off()
#---------




### System and session information
#---------
cat("\n\n\n\n... SYSTEM INFORMATION ...\n")
INFORMATION <- Sys.info()
print(as.data.frame(INFORMATION))

cat("\n\n\n\n... SESSION INFORMATION ...\n")
a<-sessionInfo()
#---------


