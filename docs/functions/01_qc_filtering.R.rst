01_qc_filtering.R
=================



### Run quality control on the raw dataset
```bash
Rscript $script_path/01_Seurat_QC.R \
    -i $main/analysis/1-QC_and_Filtering/Raw_Seurat_Object.rds \
    -c $var_to_plot \
    -s 'mouse' \
    -p $script_path/../seurat_cell_cycle \
    -f $script_path/inst_packages.R \
    -o $main/analysis/1-QC_and_Filtering \
    2>&1 | tee $main/analysis/1.QClog.txt
```
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
