The current `scRNAseq workflow` is composed of 4 main scripts, some
supplementary files for cell cycle identification. Those files are
usually located in a separate script folder and should not be altered.
The workflow is run via `run_workflow.sh` file that can be run either
from your local computer or directly on UPPMAX via `sbatch`.
 
However, the scripts are divided in parts where human choices need to be
done (such as choosing the clustering that best defines the biological
phenotype). For this purpose, one should run first run the scripts one
per time. To do so, just comment out all other scripts not to be run and
run `run_workflow.sh` via sbatch. These scripts can also be run from the
command line if necessary, but not recommended.

1.  The workflow already inludes the option to import 10X data using the
    `00_Create_Seurat_object.R` function. SMARTseq2 can be imported into
    the pipeline too, but for that you need to create a seurat object
    with the raw counts and the metadata separatelly and save as .rds
    file.

2.  Once the Seurat object is created, just run the `01_Seurat_QC.R`
    script (by commenting all other scripts and run `run_workflow.sh`
    via sbatch). It will generate several quality control plots and cell
    cycle classification scores, finaly resulting in the
    `Raw_Seurat_Object.rds` containig all the QC information, but not
    filtered. After this step, the script will also set the 1% and 99%
    thresholds for all main QC parameters and output a
    `Raw_Seurat_Object.rds`. You can choose either of these output
    objects as input for the next step. If you already have a Seurat
    object from another pipeline, it is recommended to re-run this
    script, just so all metadata names used in the next scripts are the
    same. In that case, just continue using the `Raw_Seurat_Object.rds`
    file outputed from this script, if you already have filtered cells
    previously.

3.  Once the step above is done, just run the `02_Clustering.R` script
    (by commenting all other scripts and run `run_workflow.sh` via
    sbatch). This step will compute variable genes, perform data
    scaling, PCA and tSNE for your data. After that it will screen four
    clustering methods (SNN, HDBSCAN, DBSCAN and FlowPeaks) and
    parameters and output several tSNE plots in the respective folders.
    The clustering should be inspected manually. If the tSNE coordinate
    file is already in the output folder, the tSNE will be skipped.
    Overall, they should agree with the tSNE and be in accordance with
    the differential expression in the next step.

4.  The third script will copmute differentially expressed genes
    (markers) among all clusters. In addition to this anlaysis, it also
    computes differential expression per group, instead of per cluster.
    If you have a metadata variable that you would like to compare for
    each cluster, then you can also just specify the parameter `-m` with
    the column name in your metadata. This is usefull, for example, in
    case your data has 5 clusters and you want to compare the
    differential expression between time\_points/sample\_groups within
    each cluster. Multiple arguments can be parsed comma separated, and
    will be outputed in different folders. If after the main
    differential expression you notice that 2 clusters are too similar,
    you should re-run the anlaysis using another clustering
    specification in that those clusters are unified. If differential
    expression tables are found, they will be loaded and the
    differential expression for them will be skipped. Several violin and
    tSNE plots are outputed from this step.

### Sbatch configurations

    #! /bin/bash -l
    #SBATCH -A snic2017-7-128
    #SBATCH -p core
    #SBATCH -n 2
    #SBATCH -J scRNAseq_example
    #SBATCH -t 01:00:00
    #SBATCH --mail-user myemail@nbis.se
    #SBATCH --mail-type=END

------------------------------------------------------------------------

### Load modules on UPPMAX

    module load bioinfo-tools
    module load R/3.5.0
    module load R_packages/3.5.0

This will load the latest stable version of R and its packages.

------------------------------------------------------------------------

### Define common variables and folder here
```bash
    script_path='/proj/uppstore2017171/devop/scRNAseq_pipeline/rscripts'
    main='/proj/uppstore2017171/staff/paulo/singlecell_example'
    cd $main

    var_to_plot='SampleID,batch,days_post_infection,sequencing_run'
    var_to_regress='nUMI,SampleID,percent.mito,S.Score,G2M.Score,percent.Rpl,percent.Rps,batch,sequencing_run'
```
Here is where you can define your variables that you should change:

`script_path`: the path to the scripts of the curated scRNAseq workflow.
You should not change that.

`main`: the path to your main project, where your anlaysis will be run.
In that folder, it is recommended to have your raw data in a folder
`data`. There you should also have a .csv file with the metadata for
each library. If you want to already start from a Seurat object, just
rename it and put in the folder
`$main/analysis/1-QC_and_Filtering/Raw_Seurat_Object.rds`.

`var_to_plot`: the variables from your metadata that you would like to
plot and see results separated for. E.g.: color the tSNE and barplots
based on those parameters.

`var_to_regress`: the variables from your metadata that you would like
to regress out from your data. This includes batches, embrionic time
points and other factors that might impcat on yout cell clustering. If
no batches are present, just input
`nUMI,SampleID,percent.mito,S.Score,G2M.Score`. The
`percent.Rpl,percent.Rps` are also available, and are particularly
usefull when working with cells with small RNA content, such as T cells
and ILCs. This is because these gene family impacts a lot on cell
clustering. you can experiment and see which one results in better
clustering.

------------------------------------------------------------------------

### Create Seurat object from 10x raw UMI counts

    Rscript $script_path/00_Create_Seurat_object.R \
            -i $main/data/cellranger \
            -m $main/data/metadata.csv \
            -c 'SampleID,batch,days_post_infection,sequencing_run' \
            -f $script_path/inst_packages.R \
            -o $main/analysis/1-QC_and_Filtering \
            2>&1 | tee $main/analysis/0.Import10Xlog.txt

`-i`: the input PATH with 10X files. Each sample is a folder, with the
matrix and indexes in it.

`-m`: the metadata .csv file for every 10X library, containing the
experimental or group metadata. Each sample in a row.

`-c`: the columns names from the metadata that you would like to import
in your Seurat Object. Multiple arguments are parsed comma separated.

`-f`: the path to custom scripts shared across the whole pipeline. These
are already supplied in the workflow script folder.

`-o`: the output folder. It will be created if it does not exist.

`2>&1 | tee`: the run log file.

------------------------------------------------------------------------

### Run quality control on the raw dataset

    Rscript $script_path/01_Seurat_QC.R \
            -i $main/analysis/1-QC_and_Filtering/Raw_Seurat_Object.rds \
            -c $var_to_plot \
            -s 'mouse' \
            -p $script_path/../seurat_cell_cycle \
            -f $script_path/inst_packages.R \
            -o $main/analysis/1-QC_and_Filtering \
            2>&1 | tee $main/analysis/1.QClog.txt

`-i`: the input Seurat object FILE.

`-c`: the columns names from the metadata that you would like to import
in your Seurat Object. Multiple arguments are parsed comma separated.

`-s`: the species used for single cell. This will just convert the gene
names to allow the Seurat cell cycle estimation (which uses human
symbols). Only “mouse” and “human” is supported at the moment.

`-p`: the path to Seurat cell scoring files. These are already supplied
in the workflow script folder.

`-f`: the path to custom scripts shared across the whole pipeline. These
are already supplied in the workflow script folder.

`-o`: the output folder. It will be created if it does not exist.

`2>&1 | tee`: the run log file.

------------------------------------------------------------------------

### Run clustering for the all cells

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

------------------------------------------------------------------------

### Run differential expression for the main cell types and for days\_post\_infection per cluster

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

------------------------------------------------------------------------

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

    Rscript $script_path/02_Clustering.R \
	        -i $main/analysis/2-Clustering/Clustered_Seurat_object.rds \
	        -c $var_to_plot \
    	    -r $var_to_regress \
    	    -s 'hdbscan.20,5'\
    	    -f $script_path/inst_packages.R \
    	    -o $main/analysis/3-Epithelial \
    	    2>&1 | tee $main/analysis/3.Epithelial_Clusteringlog.txt

