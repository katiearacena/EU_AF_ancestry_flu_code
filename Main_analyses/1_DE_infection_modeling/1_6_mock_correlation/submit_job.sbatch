#!/bin/bash

############################################
#SBATCH --job-name=mock_NI_corr
#SBATCH --output=logs/corr-%j.out
#SBATCH --error=logs/corr-%j.err
#SBATCH --time=03:00:00
#SBATCH --partition=lbarreiro
#SBATCH --account=pi-\lbarreiro
#SBATCH --mem=80G
############################################

module load R
Rscript script_to_corr.R

wait

echo " COMPLETED AT `date`"

## EOF ##
