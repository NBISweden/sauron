#! /bin/bash -l
this_file_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo $this_file_path



##############################
### CREATE ANALYSIS FOLDER ###
##############################
script_path=$this_file_path/'../../../scripts'
echo script_path

main=~/Downloads/sauron_tutorial_PBMC
mkdir ~/Downloads/sauron_tutorial_PBMC
cd ~/Downloads/sauron_tutorial_PBMC
echo $main

mkdir data
mkdir analysis
mkdir log

cp $this_file_path/metadata.csv data/metadata.csv



#####################
### DOWNLOAD DATA ###
#####################
# cd data
# mkdir pbmc_1k_v2
# mkdir pbmc_1k_v3
# mkdir pbmc_1k_p3
#
# curl -o pbmc_1k_v2/pbmc_1k_v2_filtered_feature_bc_matrix.h5 -O \
# http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_filtered_feature_bc_matrix.h5
# curl -o pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5 -O \
# http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5
# curl -o pbmc_1k_p3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5 -O \
# http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5
# cd ..



######################
### ACTIVATE CONDA ###
######################
source activate Sauron.v1


########################
### DEFINE VARIABLES ###
########################
var_to_plot='dataset,chemistry'
var_to_regress='nFeature_RNA,nCount_RNA,percent_mito,S.Score,G2M.Score'



#####################
### LOAD DATASETS ###
#####################
# Rscript $script_path/00_load_data.R \
# --input_path $main/'data' \
# --dataset_metadata_path $main/'data/metadata.csv' \
# --assay 'RNA' \
# --output_path $main/'analysis/1_qc' \
# 2>&1 | tee $main/log/'00_load_data_log.txt'



###########################
### RUN QUALITY CONTROL ###
###########################
# Rscript $script_path/01_qc_filter.R \
# --Seurat_object_path $main/'analysis/1_qc/raw_seurat_object.rds' \
# --columns_metadata $var_to_plot \
# --species_use 'hsapiens' \
# --remove_non_coding 'True' \
# --plot_gene_family 'RPS,RPL,mito,HB' \
# --remove_gene_family 'mito' \
# --min_gene_count '5' \
# --min_gene_per_cell '200' \
# --output_path $main/analysis/1_qc \
# 2>&1 | tee $main/log/'01_QC_log.txt'



##############################################################
### RUN DATA INTEGRATION, NORMALIZE AND GET VARIABLE GENES ###
##############################################################
# Rscript $script_path/02_integrate.R \
# --Seurat_object_path $main/'analysis/1_qc/filt_seurat_object.rds' \
# --columns_metadata $var_to_plot \
# --regress $var_to_regress \
# --var_genes 'seurat' \
# --integration_method 'mnn,dataset' \
# --cluster_use 'NONE' \
# --assay 'RNA' \
# --output_path $main/'analysis/2_clustering' \
# 2>&1 | tee $main/log/'02_integrate_log.txt'



###################################################
### RUN DIMENSIONALITY REDUCTION AND CLUSTERING ###
###################################################
Rscript $script_path/03_dr_and_cluster.R \
--Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
--columns_metadata $var_to_plot \
--regress $var_to_regress \
--PCs_use 'var,1' \
--var_genes 'seurat' \
--dim_reduct_use 'umap' \
--cluster_use 'none' \
--cluster_method 'louvain,hc,kmeans' \
--assay 'mnn' \
--output_path $main/'analysis/2_clustering' \
2>&1 | tee $main/log/'03_dr_and_cluster_log.txt'



########################################
### RUN CLUSTER CORRELATION ANALYSIS ###
########################################
# Rscript $script_path/'05_cluster_correlation.R' \
# --Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
# --clustering_use 'HC_100' \
# --exclude_cluster 'NONE' \
# --merge_cluster '0.95,0.9,0.85,0.8,0.75,0.7' \
# --output_path $main/'analysis/2_clustering/cluster_correlations' \
# 2>&1 | tee $main/'log/4_clust_corr.txt'



###################################
### RUN DIFFERENTIAL EXPRESSION ###
###################################
# Rscript $script_path/04_diff_gene_expr.R \
# --Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
# --clustering_use 'HC_12' \
# --metadata_use 'tech' \
# --exclude_cluster 'NONE' \
# --assay 'RNA' \
# --o $main/'analysis/3_diff_expr' \
# 2>&1 | tee $main/'log/4_diff_expr_log.txt'



################################################
### RUN LIGAND-RECEPTOR INTERACTION ANALYSIS ###
################################################
# Rscript $script_path/06_lig_rec_interactome.R \
# --objects_paths $main/'analysis/2_clustering/seurat_object.rds' \
# --object_names 'all_cells' \
# --object_clusters 'HC_12,1,2,3,4' \
# --lig_recp_database 'DEFAULT' \
# --ligand_objects 'all_cells' \
# --receptor_objects 'all_cells' \
# --species_use 'hsapiens' \
# --metadata_ligands 'tech' \
# --metadata_receptor 'tech' \
# --filter_thresholds '0.1,0.1,3' \
# --output_path $main/'analysis/5_Lig_Rec_interaction' \
# --assay 'RNA' \
# 2>&1 | tee $main/'log/5_Interactome_EPI_log.txt'


conda deactivate
