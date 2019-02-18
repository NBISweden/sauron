00_load_data.R
==============


### Create Seurat object from 10x raw UMI counts
```bash
Rscript $script_path/00_Create_Seurat_object.R \
    -i $main/data/cellranger \
    -m $main/data/metadata.csv \
    -c 'SampleID,batch,days_post_infection,sequencing_run' \
    -f $script_path/inst_packages.R \
    -o $main/analysis/1-QC_and_Filtering \
    2>&1 | tee $main/analysis/0.Import10Xlog.txt
```
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
