#!/bin/bash

############################################
#SBATCH --job-name=cisreg_matrix
#SBATCH --output=logs/cisreg_matrix-%j.out
#SBATCH --error=logs/cisreg_matrix-%j.err
#SBATCH --time=3:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=48G
#SBATCH --nodes=1
############################################

DATA=${1}

module load R

Rscript get_residuals.R $DATA

wait

echo " COMPLETED AT `date`"

## EOF ##