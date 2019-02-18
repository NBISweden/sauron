02_clustering.R
===============


### Run clustering for the all cells
Checkout the latest version of this repository::
    Rscript $script_path/02_Clustering.R \
        -i $main/analysis/1-QC_and_Filtering/Filt_Seurat_Object.rds \
        -c $var_to_plot \
        -r $var_to_regress \
        -s 'ClusteringName,ClusterID'\
        -f $script_path/inst_packages.R \
        -o $main/analysis/2-Clustering \
        2>&1 | tee $main/analysis/2.Clusteringlog.txt

`-i`: the input Seurat object FILE.

`-c`: the columns names from the metadata that you would like to import
in your Seurat Object. Multiple arguments are parsed comma separated.

`-r`: the columns names from the metadata that you would like to regress
out from the tSNE and clustering. Multiple arguments are parsed comma
separated.

`-s`: The clustering name and cluster ID to to use for clustering. If
either the clustering name or ID is not found, it will use the whole
data. This is usefull if you want to re-run the anlaysis on a particular
cell subset.

`-p`: the path to Seurat cell scoring files. These are already supplied
in the workflow script folder.

`-f`: the path to custom scripts shared across the whole pipeline. These
are already supplied in the workflow script folder.

`-o`: the output folder. It will be created if it does not exist.

`2>&1 | tee`: the run log file.






### What if I want to run the clustering and analysis on a specific cell population

For such cases, you can just run the `02_Clustering.R` script adjusting
the `-s` paramater the the clustering method and the respective cell
cluster to do the anlaysis. In the case below, we will run the anlaysis
on a CLustered Seurat Object. We choose the clustering method
`hdbscan.20` and the cluster `5`. Next, you should also change the
output folder `-o` and `2>&1 | tee` to something like
`$main/analysis/3-Epithelial`, so it doesn’t overithe the folder and
files from the previous step.

Please note that this is no restricted only to the selection of clusters
only, but any column in the metadata table with its respective selective
value.
```bash
Rscript $script_path/02_Clustering.R \
    -i $main/analysis/2-Clustering/Clustered_Seurat_object.rds \
    -c $var_to_plot \
    -r $var_to_regress \
    -s 'hdbscan.20,5'\
    -f $script_path/inst_packages.R \
    -o $main/analysis/3-Epithelial \
    2>&1 | tee $main/analysis/3.Epithelial_Clusteringlog.txt
```
