


library(data.table)

#PANGALO
cell_markers <- fread('https://panglaodb.se/markers/PanglaoDB_markers_25_Aug_2019.tsv.gz')


#Filter to use only conserved markers across species
cutoff <- 0.15
sel <- ((cell_markers$specificity_human^2 + cell_markers$specificity_mouse^2) ) <= cutoff^2
sel2 <- cell_markers$`official gene symbol` == "CD4"

plot(cell_markers$specificity_human, cell_markers$specificity_mouse,cex=.1,ylim=c(0,1),xlim=c(0,1),
     las=1,xlab="specificity_human",ylab="specificity_mouse",xaxs="i",yaxs="i")
points(cell_markers$specificity_human[sel], cell_markers$specificity_mouse[sel],pch=16,cex=.5,col="red")
points(cell_markers$specificity_human[sel2], cell_markers$specificity_mouse[sel2],pch=16,cex=.8,col="blue")
abline(h=.2,v=.2,lty=2,col="gray50")



#Get the markers for each cell population
filt_cell_markers <- cell_markers[sel,]
markers <- lapply( unique(cell_markers$`cell type`),function(x){
  return( casefold( cell_markers$`official gene symbol`[cell_markers$`cell type` == x] ) )
})
names(markers) <- unique(cell_markers$`cell type`)


#Save gene list to the 'cell_marker' folder in 'support_files'
rep_markers <- lapply(markers,function(x) {
  temp <- rep("",max(sapply(markers,length)))
  temp[1:length(x)] <- x
  return(temp)
} )
write.csv2(rep_markers, "~/Box/repos/sauron/support_files/cell_markers/PANGALODB_20190829.csv",row.names = F)






#CellMarker
cell_markers2 <- fread('http://biocc.hrbmu.edu.cn/CellMarker/download/all_cell_markers.txt')
















