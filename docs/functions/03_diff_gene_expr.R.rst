03_diff_gene_expr.R
===================


### Run differential expression for the main cell types::
    Rscript $script_path/03_Find_Markers.R \
        -i $main/analysis/2-Clustering/Clustered_Seurat_object.rds \
        -c 'hdbscan.17' \
        -m 'days_post_infection' \
        -e 'ExcludeCluster' \
        -f $script_path/inst_packages.R \
        -o $main/analysis/2-Clustering/DGE_per_cluster \
        2>&1 | tee $main/analysis/3-Find_Markers.txt

`-i`: the input Seurat object FILE.

`-c`: The clustering name to use for differential expression.

`-m`: If you have a metadata variable that you would like to compare for
each cluster, than you can also just put it here. Let’s say you find 5
clusters from the previous script and you want to compare the
differential expression between time\_points/sample\_groups within each
cluster. Multiple arguments are parsed comma separated.

`-e`: cluster names to exclude from the analsysis. Multiple arguments
are parsed comma separated. If the cluster is not found, it will be just
ignored. If you want to include all, just write anything that is not
your cluster.

`-f`: the path to custom scripts shared across the whole pipeline. These
are already supplied in the workflow script folder.

`-o`: the output folder. It will be created if it does not exist.

`2>&1 | tee`: the folder for the run log file.
