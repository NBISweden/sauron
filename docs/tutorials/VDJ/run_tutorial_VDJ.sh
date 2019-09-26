#! /bin/bash -l
this_file_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo $this_file_path



##############################
### CREATE ANALYSIS FOLDER ###
##############################
script_path=$this_file_path/'../../../scripts'
echo script_path

main=~/Downloads/sauron_tutorial_VDJ
mkdir ~/Downloads/sauron_tutorial_VDJ
cd ~/Downloads/sauron_tutorial_VDJ
echo $main

mkdir data
mkdir analysis
mkdir log

cp $this_file_path/metadata.csv data/metadata.csv



#####################
### DOWNLOAD DATA ###
#####################
# cd data
#
# mkdir gene_expression
# mkdir tcr_data
# mkdir ig_data
#
# #download gene expression data
# cd gene_expression
# mkdir c57; cd c57
# curl -O http://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_mm_c57bl6_pbmc_5gex/vdj_v1_mm_c57bl6_pbmc_5gex_filtered_feature_bc_matrix.h5
# cd ..
# mkdir balb; cd balb
# curl -O http://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_mm_c57bl6_pbmc_5gex/vdj_v1_mm_c57bl6_pbmc_5gex_filtered_feature_bc_matrix.h5
# cd ..
#
# #download TCR data
# cd ../tcr_data
# mkdir c57; cd c57
# curl -O http://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_mm_c57bl6_pbmc_t/vdj_v1_mm_c57bl6_pbmc_t_filtered_contig_annotations.csv
# cd ..
# mkdir balb; cd balb
# curl -O http://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_mm_balbc_pbmc_t/vdj_v1_mm_balbc_pbmc_t_filtered_contig_annotations.csv
# cd ..
#
# #download Ig data
# cd ../ig_data
# mkdir c57; cd c57
# curl -O http://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_mm_c57bl6_pbmc_b/vdj_v1_mm_c57bl6_pbmc_b_filtered_contig_annotations.csv
# cd ..
# mkdir balb; cd balb
# curl -O http://cf.10xgenomics.com/samples/cell-vdj/3.0.0/vdj_v1_mm_balbc_pbmc_b/vdj_v1_mm_balbc_pbmc_b_filtered_contig_annotations.csv
cd $main



######################
### ACTIVATE CONDA ###
######################
source activate Sauron.v1


########################
### DEFINE VARIABLES ###
########################
var_to_plot='dataset,mouse'
var_to_regress='nFeature_RNA,nCount_RNA,percent_mito,S.Score,G2M.Score'



#####################
### LOAD DATASETS ###
#####################
# Rscript $script_path/00_load_data.R \
# --input_path $main/'data/gene_expression' \
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
# --species_use 'mmusculus' \
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
# --var_genes 'scran' \
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
# --var_genes 'seurat' \
# --dim_reduct_use 'umap' \
# --cluster_use 'none' \
# --cluster_method 'louvain,hc,hdbscan' \
# --assay 'mnn' \
# --output_path $main/'analysis/2_clustering' \
# 2>&1 | tee $main/log/'03_dr_and_cluster_log.txt'



########################################
### RUN CLUSTER CORRELATION ANALYSIS ###
########################################
# Rscript $script_path/'05_cluster_correlation.R' \
# --Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
# --clustering_use 'HC_14' \
# --exclude_cluster 'NONE' \
# --merge_cluster '0.95,0.9,0.85,0.8,0.75,0.7' \
# --output_path $main/'analysis/2_clustering/cluster_correlations' \
# 2>&1 | tee $main/'log/4_clust_corr.txt'



###################################
### RUN DIFFERENTIAL EXPRESSION ###
###################################
# Rscript $script_path/04_diff_gene_expr.R \
# --Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
# --clustering_use 'HC_14' \
# --metadata_use 'dataset' \
# --exclude_cluster 'NONE' \
# --assay 'rna' \
# --output_path $main/'analysis/2_clustering/diff_expr' \
# 2>&1 | tee $main/'log/5_diff_expr_log.txt'



############################
### CELL TYPE PREDICTION ###
############################
# Rscript $script_path/cell_type_prediction.R \
# --Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
# --marker_lists $this_file_path/'../../../support_files/cell_markers/main_cell_types.csv' \
# --cluster_use 'leiden_res0.7' \
# --assay 'rna' \
# --output_path $main/'analysis/2_clustering/cell_type_prediction' \
# 2>&1 | tee $main/'log/cell_type_prediction_log.txt'



########################
### RUN VDJ ANALYSIS ### - Tcr only
########################
# Rscript $script_path/VDJ_analysis.R \
# --Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
# --VDJ_annotation_path $main/'data/tcr_data' \
# --columns_metadata 'mouse' \
# --top_TCRs '10' \
# --paired_only 'true' \
# --only_coding_cdr3 'true' \
# --same_scale 'true' \
# --assay 'rna' \
# --output_path $main/'analysis/5_vdj_analysis_tcr' \
# 2>&1 | tee $main/log/'05_vdj_analysis_tcr_log.txt'



########################
### RUN VDJ ANALYSIS ### - Ig only
########################
# Rscript $script_path/VDJ_analysis.R \
# --Seurat_object_path $main/'analysis/2_clustering/seurat_object.rds' \
# --VDJ_annotation_path $main/'data/ig_data' \
# --columns_metadata 'mouse' \
# --top_TCRs '10' \
# --paired_only 'true' \
# --only_coding_cdr3 'true' \
# --same_scale 'true' \
# --assay 'rna' \
# --output_path $main/'analysis/5_vdj_analysis_ig' \
# 2>&1 | tee $main/log/'05_vdj_analysis_ig_log.txt'



################################################
### RUN LIGAND-RECEPTOR INTERACTION ANALYSIS ###
################################################
Rscript $script_path/06_lig_rec_interactome.R \
--objects_paths $main/'analysis/2_clustering/seurat_object.rds',$main/'analysis/2_clustering/seurat_object.rds' \
--object_names 'B_cell,T_cell' \
--object_clusters 'HC_5,2;HC_5,4' \
--lig_recp_database 'DEFAULT' \
--ligand_objects 'B_cell,T_cell' \
--receptor_objects 'B_cell,T_cell' \
--species_use 'mmusculus' \
--metadata_ligands 'dataset' \
--metadata_receptor 'dataset' \
--filter_thresholds '0.1,0.1,3' \
--output_path $main/'analysis/5_Lig_Rec_interaction' \
--assay 'rna' \
2>&1 | tee $main/'log/5_Interactome_EPI_log.txt'


conda deactivate
