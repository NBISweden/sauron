#! /bin/bash -l
#SBATCH -A proj_number
#SBATCH -p core
#SBATCH -n 8
#SBATCH -J proj_ID
#SBATCH -t 16:00:00
#SBATCH --mail-user username@email.com
#SBATCH --mail-type=END


#Load modules on UPPMAX
#Uncomment the line below if working on your local computer
#module load bioinfo-tools
#module load R/3.5.0
#module load R_packages/3.5.0


#Define common variables and folder here
var_to_plot='libraryID,SampleName,SampleGroup,Batch'
var_to_regress='nUMI,percent.mito,S.Score,G2M.Score,Batch'
script_path='/proj/uppstore2017234/private/DataFromWABI/Paulo/analysis2/rscripts'
main='/proj/uppstore2017234/private/DataFromWABI/Paulo/analysis2'
cd $main



###Create Seurat object from 10x raw UMI counts
Rscript $script_path/00_load_data.R \
	-i $main/data \
	-m $main/data/metadata.csv \
	-c $var_to_plot \
	-f $script_path/inst_packages.R \
	-o $main/analysis/1-QC_and_Filtering \
	2>&1 | tee $main/analysis/0_Import10Xlog.txt



###Run quality control on all 10X samples
Rscript $script_path/01_qc_filtering.R \
	-i $main/analysis/1-QC_and_Filtering/Raw_Seurat_Object.rds \
	-c $var_to_plot \
	-s 'mouse' \
	-p $script_path/seurat_cell_cycle \
	-f $script_path/inst_packages.R \
	-o $main/analysis/1-QC_and_Filtering \
	2>&1 | tee $main/analysis/1_QClog.txt



###Run clustering for the main cell types
Rscript $script_path/02_clustering.R \
	-i $main/analysis/1-QC_and_Filtering/Filt_Seurat_Object.rds \
	-c $var_to_plot \
	-r $var_to_regress \
	-s 'ClusteringName,ClusterID'\
	-f $script_path/inst_packages.R \
	-o $main/analysis/2_Clustering \
	2>&1 | tee $main/analysis/2_Clusteringlog.txt



###Run differential expression for all main cell types
Rscript $script_path/03_diff_gene_expr.R \
	-i $main/analysis3/2-Clustering/Seurat_object.rds \
	-c 'hdbscan.19' \
	-m 'SampleGroup' \
	--exclude_cluster "0" \
	-f $script_path/inst_packages.R \
	-o $main/analysis3/3-Find_Markers \
	2>&1 | tee $main/analysis3/3.Find_Markers.txt



###Run clustering for a specific cluster (or group of clusters)
Rscript $script_path/02_clustering.R \
	-i $main/analysis3/2-Clustering/Seurat_object.rds \
	-c $var_to_plot \
	-r $var_to_regress \
	-s 'res.0.7,12'\
	-f $script_path/inst_packages.R \
	-o $main/analysis3/4-MyCluster \
	2>&1 | tee $main/analysis3/4-MyCluster_log.txt



###Run differential expression for the Fibroblast cell population
Rscript $script_path/03_diff_gene_expr.R \
	-i $main/analysis3/4-Fibroblast/Seurat_object.rds \
 	-c 'res.0.2' \
 	-m 'SampleGroup' \
 	-e '6' \
 	-f $script_path/inst_packages.R \
 	-o $main/analysis/4-MyCluster/DGE_per_cluster_res.0.2 \
 	2>&1 | tee $main/analysis3/4.Find_Markers_MyCluster.txt



###Run cluster correlation analysis
Rscript $script_path/04_cluster_correlations.R \
	-i $main/analysis/2-Clustering/Seurat_object.rds \
	-c 'flowpeaks.0.3' \
	--exclude_cluster "0" \
	-f $script_path/inst_packages.R \
	-o $main/analysis/3-Find_Markers \
	2>&1 | tee $main/analysis/3.Find_Markers.txt



###Run ligand-receptor interactome among clusters
Rscript $script_path/05_lig_rec_interactome.R \
	--objects_paths ${main}/analysis/2-Clustering/Seurat_object.rds,${main}/analysis/2-Clustering/Seurat_object.rds \
	--object_names 'Tcells,Monocytes' \
	--object_clusters 'hdbscan.99,1;hdbscan.99,2' \
	--lig_recp_database ${main}/support_files/ligand_receptor.csv \
	--ligand_objects 'Monocytes' \
	--receptor_objects 'Tcells' \
	--species_use 'human' \
	--aux_functions_path $script_path/inst_packages.R \
	--output_path $main/analysis/5-interactome \
	2>&1 | tee $main/2.Clusteringlog.txt











