05_lig_rec_interaction.R
========================



###Run ligand-receptor interactome among clusters:
.. code:: bash
    Rscript $script_path/05_lig_rec_interactome.R \
    	--objects_paths ${main}/analysis/2-Clustering/Seurat_object.rds,${main}/analysis/2-Clustering/Seurat_object.rds \
    	--object_names 'Tcells,Monocytes' \
    	--object_clusters 'hdbscan.99,1;hdbscan.99,2' \
    	--lig_recp_database ${main}/support_files/ligand_receptor.csv \
    	--ligand_objects 'Monocytes' \
    	--receptor_objects 'Tcells' \
    	--species_use 'human' \
    	--aux_functions_path $script_path/inst_packages.R \
    	--output_path $main/analysis/5-interactome \
    	2>&1 | tee $main/2.Clusteringlog.txt

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
