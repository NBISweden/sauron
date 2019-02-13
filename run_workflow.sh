#! /bin/bash -l
#SBATCH -A proj_number
#SBATCH -p core
#SBATCH -n 8
#SBATCH -J proj_ID
#SBATCH -t 16:00:00
#SBATCH --mail-user username@email.com
#SBATCH --mail-type=END



#Load modules on UPPMAX
module load bioinfo-tools
module load R/3.5.0
module load R_packages/3.5.0


#Define common variables and folder here
var_to_plot='SampleID,Embryonic_age,Sampling_date,Sex,Sample_Name'
var_to_regress='nUMI,percent.mito,S.Score,G2M.Score,Sex,Sampling_date'

script_path='/proj/uppstore2017234/private/DataFromWABI/Paulo/analysis2/rscripts'
main='/proj/uppstore2017234/private/DataFromWABI/Paulo/analysis2'
cd $main



###Create Seurat object from 10x raw UMI counts
#Rscript $script_path/00_Create_Seurat_object.R \
#	-i $main/data \
#	-m $main/data/metadata.csv \
#	-c $var_to_plot \
#	-f $script_path/inst_packages.R \
#	-o $main/analysis3/1-QC_and_Filtering \
#	2>&1 | tee $main/analysis3/0.Import10Xlog.txt



###Run quality control on all 10X samples
#Rscript $script_path/01_Seurat_QC.R \
#	-i $main/analysis3/1-QC_and_Filtering/Raw_Seurat_Object.rds \
#	-c $var_to_plot \
#	-s 'mouse' \
#	-p $script_path/seurat_cell_cycle \
#	-f $script_path/inst_packages.R \
#	-o $main/analysis3/1-QC_and_Filtering \
#	2>&1 | tee $main/analysis3/1.QClog.txt



###Run clustering for the main cell types
#Rscript $script_path/02_Clustering.R \
#	-i $main/analysis3/1-QC_and_Filtering/Filt_Seurat_Object.rds \
#	-c $var_to_plot \
#	-r $var_to_regress \
#	-s 'ClusteringName,ClusterID'\
#	-f $script_path/inst_packages.R \
#	-o $main/analysis3/2-Clustering \
#	2>&1 | tee $main/analysis3/2.Clusteringlog.txt



###Run differential expression for the main cell types and
###for days_post_infection per cluster
Rscript $script_path/03_Find_Markers.R \
		-i $main/analysis3/2-Clustering/Seurat_object.rds \
		-c 'hdbscan.19' \
		-m 'Embryonic_age' \
		--exclude_cluster "0" \
		-f $script_path/inst_packages.R \
		-o $main/analysis3/3-Find_Markers \
		2>&1 | tee $main/analysis3/3.Find_Markers.txt



###Run clustering for the Fibroblasts
Rscript $script_path/02_Clustering.R \
	-i $main/analysis3/2-Clustering/Seurat_object.rds \
	-c $var_to_plot \
	-r $var_to_regress \
	-s 'hdbscan.19,12'\
	-f $script_path/inst_packages.R \
	-o $main/analysis3/4-Fibroblast \
	2>&1 | tee $main/analysis3/4-Fibroblast_log.txt



###Run differential expression for the Fibroblast cell population
#Rscript $script_path/03_Find_Markers.R \
#	-i $main/analysis3/4-Fibroblast/Seurat_object.rds \
# 	-c 'res.0.2' \
# 	-m 'Embryonic_age' \
# 	-e '6' \
# 	-f $script_path/inst_packages.R \
# 	-o $main/analysis/4-Fibroblast/DGE_per_cluster_res.0.2 \
# 	2>&1 | tee $main/analysis3/4.Find_Markers_Fibroblast.txt


###Run clustering for the Epithelial cells
Rscript $script_path/02_Clustering.R \
	-i $main/analysis3/2-Clustering/Seurat_object.rds \
 	-c $var_to_plot \
	-r $var_to_regress \
	-s 'hdbscan.19,3,7'\
	-f $script_path/inst_packages.R \
	-o $main/analysis3/5-Epithelial \
	2>&1 | tee $main/analysis3/5.Epithelial_log.txt



#Rscript $script_path/03_Find_Markers.R \
#	-i $main/analysis/5-Epithelial/Seurat_object.rds \
#	-c 'res.0.2' \
#	-m 'Embryonic_age' \
#	-e '5' \
#	-f $script_path/inst_packages.R \
#	-o $main/analysis/5-Epithelial/DGE_per_cluster_res0.2 \
#	2>&1 | tee $main/analysis/5.Find_Markers_Epithelial_res0.2.txt


