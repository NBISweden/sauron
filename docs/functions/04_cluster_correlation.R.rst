04_cluster_correlation.R
========================

###Run cluster correlation analysis
Rscript $script_path/04_cluster_correlations.R \
	-i $main/analysis/2-Clustering/Seurat_object.rds \
	-c 'flowpeaks.0.3' \
	--exclude_cluster "0" \
	-f $script_path/inst_packages.R \
	-o $main/analysis/3-Find_Markers \
	2>&1 | tee $main/analysis/3.Find_Markers.txt
