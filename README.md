# Single cell Analysis workflow
[![Documentation Status](https://readthedocs.org/projects/single-cell-analysis/badge/?version=latest)](https://single-cell-analysis.readthedocs.io/en/latest/?badge=latest)

Paulo Czarnewski

This is a repo for easily running single cell data analysis on UPPMAX via SLURM queueing system.

Please check the complete manual for further reading:
https://single-cell-analysis.readthedocs.io/en/latest/index.html

------------------------------------------------------------------------


## SQS - Super Quick Start

The workflow consists of 3 main steps as follows:

1.First, you will need to clone this repo into your UPPMAX account.
```bash
git clone https://czarnewski@bitbucket.org/scilifelab-lts/single_cell_analysis.git
```

2.Then, just edit the `run_workflow.sh` file to match your file directories and refer to the folder containing the scripts from the repo. Please read below on how to edit the `run_workflow.sh` file.
```bash
gedit run_workflow.sh &
```

3.After all parameters and file paths have been set, just submit the job to SLURM.
```bash
sbatch run_workflow.sh
```
Once the job is submitted, just relax and take a glass of fresh coconut water.

------------------------------------------------------------------------


## Detailed documentation on the workflow

The current `scRNAseq workflow` is composed of 4 main scripts, some
supplementary files for cell cycle identification. Those files are
usually located in a separate script folder and should not be altered.
The workflow is run via `run_workflow.sh` file that can be run either
from your local computer or directly on UPPMAX via `sbatch`.

However, the scripts are divided in parts where human choices need to be
done (such as choosing the clustering that best defines the biological
phenotype). For this purpose, one should run first run the scripts one
per time. To do so, just comment out all other scripts not to be run and
run `run_workflow.sh` via `sbatch`. These scripts can also be run from the
command line if necessary, but not recommended.

1.  The workflow already includes the option to import 10X data using the
    `00_Create_Seurat_object.R` function. SMARTseq2 can be imported into
    the pipeline too, but for that you need to create a Seurat object
    with the raw counts and the metadata separately and save as `.rds`
    file.

2.  Once the Seurat object is created, just run the `01_Seurat_QC.R`
    script (by commenting all other scripts and run `run_workflow.sh`
    via `sbatch`). It will generate several quality control plots and cell
    cycle classification scores, finally resulting in the
    `Raw_Seurat_Object.rds` containing all the QC information, but not
    filtered. After this step, the script will also set the 1% and 99%
    thresholds for all main QC parameters and output a
    `Raw_Seurat_Object.rds`. You can choose either of these output
    objects as input for the next step. If you already have a Seurat
    object from another pipeline, it is recommended to re-run this
    script, just so all metadata names used in the next scripts are the
    same. In that case, just continue using the `Raw_Seurat_Object.rds`
    file outputted from this script, if you already have filtered cells
    previously.

3.  Once the step above is done, just run the `02_Clustering.R` script
    (by commenting all other scripts and run `run_workflow.sh` via
    `sbatch`). This step will compute variable genes, perform data
    scaling, PCA and tSNE for your data. After that it will screen four
    clustering methods (SNN, HDBSCAN, DBSCAN and FlowPeaks) and
    parameters and output several tSNE plots in the respective folders.
    The clustering should be inspected manually. If the tSNE coordinate
    file is already in the output folder, the tSNE will be skipped.
    Overall, they should agree with the tSNE and be in accordance with
    the differential expression in the next step.

4.  The third script will compute differentially expressed genes
    (markers) among all clusters. In addition to this analysis, it also
    computes differential expression per group, instead of per cluster.
    If you have a metadata variable that you would like to compare for
    each cluster, then you can also just specify the parameter `-m` with
    the column name in your metadata. This is useful, for example, in
    case your data has 5 clusters and you want to compare the
    differential expression between time\_points/sample\_groups within
    each cluster. Multiple arguments can be parsed comma separated, and
    will be outputted in different folders. If after the main
    differential expression you notice that 2 clusters are too similar,
    you should re-run the analysis using another clustering
    specification in that those clusters are unified. If differential
    expression tables are found, they will be loaded and the
    differential expression for them will be skipped. Several violin and
    tSNE plots are outputted from this step.

### Sbatch configurations
```bash
#! /bin/bash -l
#SBATCH -A #UPPMAX_PROJECT_ID_NUMBER
#SBATCH -p core
#SBATCH -n 2
#SBATCH -J scRNAseq_example
#SBATCH -t 01:00:00
#SBATCH --mail-user myemail@nbis.se
#SBATCH --mail-type=END
```
------------------------------------------------------------------------

### Install dependencies on UPPMAX

You can install all the required dependencies by using the workflow- and
package-managing software `Conda`, which will create a separate and contained
environment within which you can run the pipeline. This greatly increases the
reproducibility of the workflow, as exact package versions can be specified.

Start by loading Conda and creating the environment and the dependencies.

```bash
module load conda
conda env create --prefix conda_env --file environment.yml
```

Activate the environment so that all executables (`R`, `Rscript`, *etc.*) are
using the ones specified by Conda. The environment should always be active when
you are running the pipeline.

```bash
conda activate conda_env/
```

You can alternatively load the latest stable versions of each dependency
that is available on Uppmax, but this will not guarantee reproducibility
or compatibility with future changes to the packages themselves.

```bash
module load bioinfo-tools
module load R/3.5.0
module load R_packages/3.5.0
```

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
