===============
Guided Analysis
===============
.. highlight:: sh

Step 1 - Pre-requisites
-----------------------

You should first install:
- R
- git

If you have windows, you should also install:
- MobaXterm (to access unix command line)
- 7zip (to unzip tar files)


Step 1 - Setup
--------------
open Terminal or MobaXterm

create your project folder (if you don't already have one)::
mkdir ~/Desktop/MyProject

.. code-block:: bash
    $ mkdir ~/Desktop/MyProject

go to your project folder::
cd ~/Desktop/MyProject

.. code-block:: bash
    $ cd ~/Desktop/MyProject

clone the repository for your single cell analysis::
git clone https://czarnewski@bitbucket.org/scilifelab-lts/single_cell_analysis.git

Looking at directory, you should see a folder named single_cell_analysis containing the whole repository:
.. image:: ./img/tutorial/repo_folder.png
    :width: 800
    :alt: Weighted hits plot

download the data to be used in this tutorial. We will be using data from 10X Genomics::
https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

You should extract the folder contents into the 'data' directory, renaming the folder for each sample number. The final folder directory should look like this:
.. image:: ./img/tutorial/data_folder.png
    :width: 800
    :alt: Weighted hits plot


Step 2
------

Step 3
------

Step 4
------

Step 5
------
