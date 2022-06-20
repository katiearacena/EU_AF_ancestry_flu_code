#!/bin/bash

############################################
#SBATCH --job-name=WGBS_popDE
#SBATCH --output=logs/WGBS-PopDE-%j.out
#SBATCH --error=logs/WGBS-PopDE-%j.err
#SBATCH --time=36:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=256G
#SBATCH --nodes=1
############################################

# make logs directory.
mkdir logs

module load R
Rscript WGBS_PopDE.R

wait

echo " COMPLETED AT `date`"
