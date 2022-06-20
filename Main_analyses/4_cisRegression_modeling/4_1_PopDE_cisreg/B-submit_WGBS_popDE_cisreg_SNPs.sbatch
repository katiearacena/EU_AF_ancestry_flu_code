#!/bin/bash

############################################
#SBATCH --job-name=WGBS_cisregPopDE
#SBATCH --output=logs/WGBS_cisregPopDE-%j.out
#SBATCH --error=logs/WGBS_cisregPopDE-%j.err
#SBATCH --time=96:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=200G
#SBATCH --nodes=1
############################################

# make logs directory.
mkdir logs

module load R
Rscript B-WGBS_PopDE_cisreg_SNPs.R

wait

echo " COMPLETED AT `date`"
