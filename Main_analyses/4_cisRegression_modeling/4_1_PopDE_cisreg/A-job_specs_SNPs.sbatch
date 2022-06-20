#!/bin/bash

############################################
#SBATCH --job-name=cisregPopDE
#SBATCH --output=logs/cisregPopDE-%j.out
#SBATCH --error=logs/cisregPopDE-%j.err
#SBATCH --time=168:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=48G
#SBATCH --nodes=1
############################################

DATA=${1}

module load R

Rscript A-popDE_cisreg_SNPs.R $DATA

wait

echo " COMPLETED AT `date`"

## EOF ##
