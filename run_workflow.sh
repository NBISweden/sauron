#! /bin/bash -l
#SBATCH -A proj_number
#SBATCH -p core
#SBATCH -n 8
#SBATCH -J proj_ID
#SBATCH -t 16:00:00
#SBATCH --mail-user username@email.com
#SBATCH --mail-type=END



###################################
### LOAD MODULE ON HPC - UPPMAX ###
###################################
# module load bioinfo-tools
# module load R/3.5.0
# module load R_packages/3.5.0
# module load conda



##################################
### ACTIVATE CONDA ENVIRONMENT ###
##################################
# conda env create -f environment.yml
source activate Sauron.v1



########################
### DEFINE VARIABLES ###
########################
var_to_plot='Sequencing_ID,Sample_Name,Sample_ID,Batch,Group'
var_to_regress='nUMI,percent.mito,S.Score,G2M.Score'
script_path='/Users/paulo.barenco/Box/repos/single_cell_analysis/scripts'
main='/Users/paulo.barenco/Desktop/Desktop_stuff/MyProject/single_cell_analysis'
cd $main



#####################
### LOAD DATASETS ###
#####################
# Rscript $script_path/00_load_data.R \
#   --input_path $main/'data' \
#   --dataset_metadata_path $main/'data/metadata.csv' \
#   --columns_metadata $var_to_plot \
#   --integrate 'TRUE' \
#   --output_path $main/'analysis/1_qc' \
#   2>&1 | tee $main/'00_load_data_log.txt'



#######################
### QUALITY CONTROL ###
#######################
Rscript $script_path/01_qc_filter.R \
	--Seurat_object_path $main/analysis/1_qc/Raw_Seurat_Object.rds \
	--columns_metadata $var_to_plot \
	--species_use 'hsapiens' \
	--remove_non_coding 'True' \
  --plot_gene_family 'RPS,RPL,MT-,HB[AB]' \
	--remove_gene_family 'MT-' \
	--output_path $main/analysis/1_qc \
	2>&1 | tee $main/'01_QC_log.txt'



###########################
### DATASET INTEGRATION ###
###########################
# Rscript $script_path/02_integrate.R \
# 	-i $main/analysis/1_qc/Raw_Seurat_Object.rds \
# 	-c $var_to_plot \
# 	-r $var_to_regress \
# 	-p 'top,10' \
# 	-v 'cca' \
# 	-s 'CLUSTERING_NAME,CLUSTER_ID' \
# 	-m 'hdbscan,flowpeaks,snn' \
# 	-o $main/analysis/2_clustering \
# 	2>&1 | tee $main/analysis/02_clustering_log.txt



###Run clustering for the main cell types
Rscript $script_path/02_clustering.R \
    -i $main/analysis/1-QC_and_Filtering/Filt_Seurat_Object.rds \
    -c $var_to_plot \
    -r $var_to_regress \
    -m 'hdbscan' \
    -s 'ClusteringName,ClusterID'\
    -o $main/analysis/2_Clustering \
    2>&1 | tee $main/analysis/2_Clusteringlog.txt



###Run differential expression for all main cell types
Rscript $script_path/03_diff_gene_expr.R \
    -i $main/analysis/2_Clustering/Seurat_object.rds \
    -c 'hdbscan.19' \
    -m 'SampleGroup' \
    --exclude_cluster "0" \
    -o $main/analysis/3_Find_Markers \
    2>&1 | tee $main/analysis/3.Find_Markers.txt



###Run clustering for a specific cluster (or group of clusters)
Rscript $script_path/02_clustering.R \
    -i $main/analysis/2-Clustering/Seurat_object.rds \
    -c $var_to_plot \
    -r $var_to_regress \
    -s 'res.0.7,12'\
    -f $script_path/inst_packages.R \
    -o $main/analysis/4-MyCluster \
    2>&1 | tee $main/analysis/4-MyCluster_log.txt



###Run differential expression for the Fibroblast cell population
Rscript $script_path/03_diff_gene_expr.R \
    -i $main/analysis/4-Fibroblast/Seurat_object.rds \
     -c 'res.0.2' \
     -m 'SampleGroup' \
     -e '6' \
     -f $script_path/inst_packages.R \
     -o $main/analysis/4-MyCluster/DGE_per_cluster_res.0.2 \
     2>&1 | tee $main/analysis/4.Find_Markers_MyCluster.txt



###Run cluster correlation analysis
Rscript $script_path/04_cluster_correlation.R \
    -i $main/analysis/2_Clustering/Seurat_object.rds \
    -c 'hdbscan.19' \
    --exclude_cluster "0" \
    -o $main/analysis/3_Find_Markers \
    2>&1 | tee $main/analysis/3.Find_Markers.txt



###Run ligand-receptor interactome among clusters
Rscript $script_path/05_lig_rec_interactome.R \
    --objects_paths ${main}/analysis/2_Clustering/Seurat_object.rds,${main}/analysis/2_Clustering/Seurat_object.rds \
    --object_names 'Tcells,Monocytes' \
    --object_clusters 'hdbscan.19,1;hdbscan.19,2' \
    --lig_recp_database ${main}/support_files/ligand_receptor/ligand_receptor_pairs.csv \
    --ligand_objects 'Monocytes' \
    --receptor_objects 'Tcells' \
    --species_use 'human' \
    --output_path $main/analysis/5_interactome \
    2>&1 | tee $main/2.Clusteringlog.txt











