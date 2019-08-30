#! /bin/bash -l

##################################
### ACTIVATE CONDA ENVIRONMENT ###
##################################
if [[ $(conda env list) == *Sauron.v1* ]]
then
    echo 'Sauron.v1 environment was found and will be used'
else
    echo 'Sauron.v1 environment was NOT found and will be created now'
    export CONDA_ENVS_PATH=$main/Conda_env_Sauron.v1
    conda env create -n Sauron.v1 -f sauron/environment.yml
fi
source activate $main/Conda_env_Sauron.v1/Sauron.v1



############################
### DEFINE WORKFLOW PATHS ###
############################
cd $main
script_path=$main/scripts
[ ! -d analysis ] && mkdir -p analysis
[ ! -d log ] && mkdir -p log

echo "Sauron pipeline is running ..."
echo "MAIN directory is: "$main
echo "SCRIPT directory is: "$script_path
