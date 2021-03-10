#!/bin/bash

############################################
#SBATCH --job-name=mash_QTL
#SBATCH --output=logs/mash_QTL-%j.out
#SBATCH --error=logs/mash_QTL-%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=50G
#SBATCH --nodes=1
############################################

DATA=${1}

module load R/3.6.3-no-openblas

Rscript STR_mash_matrixeqtl.R $DATA

wait

echo " COMPLETED AT `date`"

## EOF ##
