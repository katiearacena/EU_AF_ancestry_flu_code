#!/bin/bash

############################################
#SBATCH --job-name=WGBS_popDR
#SBATCH --output=logs/WGBS-PopDR-%j.out
#SBATCH --error=logs/WGBS-PopDR-%j.err
#SBATCH --time=03:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=256G
#SBATCH --nodes=1
############################################

# make logs directory.
mkdir logs

module load R

Rscript WGBS_PopDR.R

wait

echo " COMPLETED AT `date`"
