# Single cell Analysis workflow
[![Documentation Status](https://readthedocs.org/projects/single-cell-analysis/badge/?version=latest)](https://single-cell-analysis.readthedocs.io/en/latest/?badge=latest)

Paulo Czarnewski

This is a repo for easily running single cell data analysis locally or via SLURM queueing system.
------------------------------------------------------------------------

## Quick Start

The workflow consists of 4 main steps as follows:

1.First, you will need to clone this repo into your project folder.
```bash
git clone https://czarnewski@bitbucket.org/scilifelab-lts/sauron.git
```

2.Add your data and metadata into the `data` directory. One dataset matrix per folder (i.e. one plate per folder or one 10X lane per folder). Name each folder as the desired sample names. The sample names should match the names in the 1st column of you metadata csv file.

3.Then, just edit the `run_workflow.sh` file to configure how you would like to run the workflow. Each script has a set of variables that you can define to run your analysis. Many of the steps have sensible thresholds and cutoff as defaults, which can be modified if needed.
```bash
gedit run_workflow.sh &
```

4.After all parameters and file paths have been set, just submit the job to SLURM.
```bash
sbatch run_workflow.sh
```
