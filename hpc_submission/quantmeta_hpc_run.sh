#!/bin/bash

###############################
#                             #
#  1) Job Submission Options  #
#                             #
###############################

# Name
#SBATCH --job-name=quantmeta

# Resources
# For MPI, increase ntasks-per-node
# For multithreading, increase cpus-per-task
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000mb
#SBATCH --time=02:00:00

# Account
#SBATCH --account=duhaimem1 ###Update to hpc specifications
#SBATCH --partition=standard ###Update to hpc specifications

# Logs
#SBATCH --mail-user=klangenf@umich.edu ###Update to be use specific
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=Logs/%x-%j.out

# Environment
##SBATCH --export=ALL

###Before running, have snakemake conda environment active by running "conda activate snakemake"
source /etc/profile.d/http_proxy.sh

# List compute nodes allocated to the job
if [[ $SLURM_JOB_NODELIST ]] ; then
    echo "Running on"
    scontrol show hostnames $SLURM_JOB_NODELIST
    echo -e "\n"
fi



#####################
#                   #
#  2) Job Commands  #
#                   #
#####################

# Initiating snakemake and running workflow in cluster mode
DIR='/home/klangenf/anaconda/envs' ### Update to user specific anaconda environment locations
snakemake --profile Config --latency-wait 300 --use-conda --conda-prefix $DIR --conda-frontend mamba --snakefile Snakefile-quantmeta
