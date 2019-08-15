#! /bin/bash -l
#SBATCH -A snic2019-3-325
#SBATCH -p core
#SBATCH -n 8
#SBATCH -J M_kasper_seurat_v3
#SBATCH -t 10:00:00
#SBATCH --mail-user paulo.czarnewski@nbis.com
#SBATCH --mail-type=END



###################################
### LOAD MODULE ON HPC - UPPMAX ###
###################################
# module load bioinfo-tools
# module load R/3.5.0
# module load R_packages/3.5.0
# module load conda


########################
### DEFINE VARIABLES ###
########################
var_to_plot='sample_id,embryonic_age,genotype,dataset,sample_no,sample_code,sex'
var_to_regress='nFeature_RNA,percent_mito,S.Score,G2M.Score'
script_path='/crex/proj/uppstore2019086/private/DataFromWABI/Paulo/sauron/scripts'
main='/crex/proj/uppstore2019086/private/DataFromWABI/Paulo'
cd $main
mkdir analysis
mkdir log


##################################
### ACTIVATE CONDA ENVIRONMENT ###
##################################
if [[ $(conda env list) == *Sauron.v1* ]]
then
    echo 'Sauron.v1 environment was found and will be used'
else
    echo 'Sauron.v1 environment was NOT found and will be created now'
    export CONDA_ENVS_PATH=$main/Conda_env_Sauron.v1
    conda env create -n Sauron.v1 -f sauron/environment.yml
fi
source activate $main/Conda_env_Sauron.v1/Sauron.v1


#####################
### LOAD DATASETS ###
#####################
 Rscript $script_path/00_load_data.R \
   --input_path $main/'data' \
   --dataset_metadata_path $main/'data/metadata.csv' \
   --columns_metadata $var_to_plot \
   --integrate 'TRUE' \
   --output_path $main/'analysis/1_qc' \
   2>&1 | tee $main/log/'00_load_data_log.txt'



###########################
### RUN QUALITY CONTROL ###
###########################
 Rscript $script_path/01_qc_filter.R \
     --Seurat_object_path $main/'analysis/1_qc/Raw_Seurat_Object.rds' \
     --columns_metadata $var_to_plot \
     --species_use 'mmusculus' \
     --remove_non_coding 'True' \
     --plot_gene_family 'RPS,RPL,mito,HB' \
     --remove_gene_family 'mito' \
     --min_gene_count '5' \
     --min_gene_per_cell '200' \
     --output_path $main/analysis/1_qc \
     2>&1 | tee $main/log/'01_QC_log.txt'



##############################################################
### RUN DATA INTEGRATION, NORMALIZE AND GET VARIABLE GENES ###
##############################################################
 Rscript $script_path/02_integrate.R \
     --Seurat_object_path $main/'analysis/1_qc/Filt_Seurat_Object.rds' \
     --columns_metadata $var_to_plot \
     --regress $var_to_regress \
     --var_genes 'seurat' \
     --integration_method 'mnn,sample_id' \
     --cluster_use 'NONE' \
     --assay 'RNA' \
     --output_path $main/'analysis/2_clustering' \
     2>&1 | tee $main/log/'02_integrate_log.txt'



####################################
### RUN DIMENSIONALITY REDUCTION ###
####################################
 Rscript $script_path/03_dr_and_cluster.R \
     --Seurat_object_path $main/'analysis/2_clustering/Seurat_Object.rds' \
     --columns_metadata $var_to_plot \
     --regress $var_to_regress \
     --PCs_use 'var,1' \
     --var_genes 'seurat' \
     --dim_reduct_use 'umap' \
     --cluster_use 'none' \
     --cluster_method 'HC,kmeans' \
     --output_path $main/'analysis/2_clustering' \
     2>&1 | tee $main/log/'03_dr_and_cluster_log.txt'



########################################
### RUN CLUSTER CORRELATION ANALYSIS ###
########################################
Rscript $script_path/'05_cluster_correlation.R' \
	--Seurat_object_path $main/'analysis/2_clustering/Seurat_object.rds' \
	--clustering_use 'HC_100' \
	--exclude_cluster 'NONE' \
	--merge_cluster '0.95,0.9,0.85,0.8,0.75,0.7' \
	--output_path $main/'analysis/2_clustering/cluster_correlations' \
	2>&1 | tee $main/'log/4_clust_corr.txt'



###################################
### RUN DIFFERENTIAL EXPRESSION ###
###################################
# Rscript $script_path/04_diff_gene_expr.R \
# 	--Seurat_object_path $main/'analysis/2_clustering/Seurat_object.rds' \
# 	--clustering_use 'HC_12' \
# 	--metadata_use 'tech' \
# 	--exclude_cluster 'NONE' \
# 	--assay 'RNA' \
# 	--o $main/'analysis/3_diff_expr' \
# 	2>&1 | tee $main/'log/4_diff_expr_log.txt'



################################################
### RUN LIGAND-RECEPTOR INTERACTION ANALYSIS ###
################################################
# Rscript $script_path/06_lig_rec_interactome.R \
# 	--objects_paths $main/'analysis/2_clustering/Seurat_object.rds' \
# 	--object_names 'all_cells' \
# 	--object_clusters 'HC_12,1,2,3,4' \
# 	--lig_recp_database 'DEFAULT' \
# 	--ligand_objects 'all_cells' \
# 	--receptor_objects 'all_cells' \
# 	--species_use 'hsapiens' \
# 	--metadata_ligands 'tech' \
# 	--metadata_receptor 'tech' \
# 	--filter_thresholds '0.1,0.1,3' \
# 	--output_path $main/'analysis/5_Lig_Rec_interaction' \
# 	--assay 'RNA' \
# 	2>&1 | tee $main/'log/5_Interactome_EPI_log.txt'


conda deactivate
