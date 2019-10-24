#! /bin/bash -l

this_file_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo $this_file_path


##############################
### CREATE ANALYSIS FOLDER ###
##############################
script_path=$this_file_path/'../../../scripts'
echo script_path

main=~/Downloads/sauron_tutorial_PBMC
mkdir $main
cd $main
echo $main

# mkdir data
# mkdir analysis
# mkdir log

# cp $this_file_path/metadata.csv data/metadata.csv



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
var_to_regress='S.Score,G2M.Score'



#####################
### LOAD DATASETS ###
#####################
# Rscript $script_path/00_load_data.R \
# --input_path $main/'data' \
# --dataset_metadata_path $main/'data/metadata.csv' \
# --assay 'rna' \
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
# --assay 'rna' \
# --output_path $main/analysis/1_qc \
# 2>&1 | tee $main/log/'01_QC_log.txt'



##############################################################
### RUN DATA INTEGRATION, NORMALIZE AND GET VARIABLE GENES ###
##############################################################
# Rscript $script_path/02_integrate.R \
# --Seurat_object_path $main/'analysis/1_qc/filt_seurat_object.rds' \
# --columns_metadata $var_to_plot \
# --regress $var_to_regress \
# --var_genes 'scran,0.001' \
# --integration_method 'mnn,dataset' \
# --cluster_use 'NONE' \
# --assay 'rna' \
# --output_path $main/'analysis/2_clustering' \
# 2>&1 | tee $main/log/'02_integrate_log.txt'



###################################################
### RUN DIMENSIONALITY REDUCTION AND CLUSTERING ###
###################################################
# Rscript $script_path/03_dr_and_cluster.R \
# --Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
# --columns_metadata $var_to_plot \
# --regress $var_to_regress \
# --PCs_use 'top,30' \
# --var_genes 'scran,0.001' \
# --dim_reduct_use 'umap' \
# --cluster_use 'none' \
# --cluster_method 'louvain,hc' \
# --assay 'mnn' \
# --output_path $main/'analysis/2_clustering' \
# 2>&1 | tee $main/log/'03_dr_and_cluster_log.txt'



########################################
### RUN CLUSTER CORRELATION ANALYSIS ###
########################################
# Rscript $script_path/'05_cluster_correlation.R' \
# --Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
# --clustering_use 'HC_17' \
# --exclude_cluster 'NONE' \
# --merge_cluster '0.95,0.9,0.85,0.8,0.75,0.7' \
# --output_path $main/'analysis/2_clustering/cluster_correlations' \
# 2>&1 | tee $main/'log/4_clust_corr.txt'



###################################
### RUN DIFFERENTIAL EXPRESSION ###
###################################
# Rscript $script_path/04_diff_gene_expr.R \
# --Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
# --clustering_use 'merged_0.85' \
# --metadata_use 'dataset' \
# --exclude_cluster 'NONE' \
# --assay 'rna' \
# --output_path $main/'analysis/2_clustering/diff_expr' \
# 2>&1 | tee $main/'log/5_diff_expr_log.txt'



############################
### CELL TYPE PREDICTION ###
############################
Rscript $script_path/cell_type_prediction.R \
--Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
--marker_lists $this_file_path/'../../../support_files/cell_markers/main_cell_types.csv' \
--cluster_use 'merged_0.95' \
--assay 'rna' \
--output_path $main/'analysis/2_clustering/cell_type_prediction' \
2>&1 | tee $main/'log/cell_type_prediction_log.txt'



##############################
### PLOT GENES OF INTEREST ###
##############################
# cp $script_path/'../support_files/cell_markers/T_cell.csv' $main/'T_cell.csv'
#
# Rscript $script_path/plot_gene_list.R \
# --Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
# --gene_list $main/'T_cell.csv' \
# --match_type 'exact' \
# --clustering_use 'merged_0.95' \
# --assay 'rna' \
# --output_path $main/'analysis/2_clustering/gene_plots' \
# 2>&1 | tee $main/'log/5_gene_plots_log.txt'



################################################
### RUN LIGAND-RECEPTOR INTERACTION ANALYSIS ###
################################################
# Rscript $script_path/06_lig_rec_interactome.R \
# --objects_paths $main/'analysis/2_clustering/seurat_object.rds' \
# --object_names 'all' \
# --object_clusters 'merged_0.75,1,2,4,7,9,16,17' \
# --lig_recp_database 'DEFAULT' \
# --ligand_objects 'all' \
# --receptor_objects 'all' \
# --species_use 'hsapiens' \
# --metadata_ligands 'none' \
# --metadata_receptor 'none' \
# --filter_thresholds '0,1,3' \
# --output_path $main/'analysis/5_Lig_Rec_interaction' \
# --assay 'rna' \
# 2>&1 | tee $main/'log/5_Interactome_EPI_log.txt'



##############################################################
### RUN DATA INTEGRATION, NORMALIZE AND GET VARIABLE GENES ### - T cells
##############################################################
# Rscript $script_path/02_integrate.R \
# --Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
# --columns_metadata $var_to_plot \
# --regress $var_to_regress \
# --var_genes 'scran,0.005' \
# --integration_method 'mnn,dataset,20' \
# --cluster_use 'merged_0.75,2' \
# --assay 'rna' \
# --output_path $main/'analysis/3_T_cells' \
# 2>&1 | tee $main/log/'02_integrate_Tcells_log.txt'



###################################################
### RUN DIMENSIONALITY REDUCTION AND CLUSTERING ### - T cells
###################################################
# Rscript $script_path/07_trajectory.R \
# --Seurat_object_path $main/'analysis/3_T_cells/seurat_object.rds' \
# --metadata_use 'louvain_0.9' \
# --reduction_use 'umap10' \
# --reduction_visualize 'umap' \
# --method_use 'ddrtree' \
# --cluster_use 'none' \
# --no_traj_components '3' \
# --no_of_paths 'none' \
# --start_cluster 'none' \
# --diff_testing 'no' \
# --assay 'rna' \
# --output_path $main/'analysis/integrated_1000_k10/clustering/traj_ddrtree' \
# 2>&1 | tee $main/'log/07_trajectory_log.txt'
#

##############################
### PLOT GENES OF INTEREST ###
##############################
# Rscript $script_path/plot_gene_list.R \
# --Seurat_object_path $main/'analysis/3_T_cells/seurat_object.rds' \
# --gene_list $main/'T_cell.csv' \
# --match_type 'exact' \
# --clustering_use 'HC_6' \
# --assay 'rna' \
# --output_path $main/'analysis/3_T_cells/gene_plots' \
# 2>&1 | tee $main/'log/5_gene_plots_log.txt'



######################
### RUN TRAJECTORY ### - T cells
######################
# Rscript $script_path/07_trajectory.R \
# --Seurat_object_path $main/'analysis/3_T_cells/seurat_object.rds' \
# --metadata_use 'HC_6' \
# --reduction_use 'umap10' \
# --reduction_visualize 'umap' \
# --method_use 'ddrtree' \
# --cluster_use 'none' \
# --no_traj_components '3' \
# --no_of_paths 'none' \
# --start_cluster '2' \
# --diff_testing 'no' \
# --assay 'rna' \
# --output_path $main/'analysis/3_T_cells/traj_ddrtree' \
# 2>&1 | tee $main/'log/07_trajectory_log.txt'











conda deactivate
