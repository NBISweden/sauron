Single cell analysis workflow
=============================

Welcome to the Single cell Analysis workflow.
The source code of this workflow can be found on `bitbucket <https://bitbucket.org/czarnewski/single_cell_analysis/src/master/>`_.

If this is your first time using the workflow, please try running the guided analysis first to have a feeling on how to use this workflow effectively.


.. toctree::
   :maxdepth: 1
   :caption: Using the workflow

   workflow/workflow
   workflow/guided_analysis



.. toctree::
  :maxdepth: 1
  :caption: Main functions

  functions/00_load_data.R
  functions/01_qc_filtering.R
  functions/02_clustering.R
  functions/03_diff_gene_expr.R
  functions/04_cluster_correlation.R


.. toctree::
  :maxdepth: 1
  :caption: Support functions

  functions/inst_packages.R
  functions/plot_gene_list.R



.. toctree::
  :maxdepth: 1
  :caption: Misc

  source/license
  source/contact
  source/citing
