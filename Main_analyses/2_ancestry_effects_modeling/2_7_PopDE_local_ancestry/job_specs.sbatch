#!/bin/bash

############################################
#SBATCH --job-name=LA_popDE_jobs
#SBATCH --output=logs/LA_PopDE-%j.out
#SBATCH --error=logs/LA_PopDE-%j.err
#SBATCH --account=pi-\lbarreiro
#SBATCH --time=128:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=24G
#SBATCH --nodes=1
############################################

DATA=${1}

module load gsl
module load R

Rscript PopDE_LA.R $DATA FALSE

wait

echo " COMPLETED AT `date`"

## EOF ##
