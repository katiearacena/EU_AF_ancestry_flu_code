#!/bin/bash

############################################
#SBATCH --job-name=SNP_STR_PVE
#SBATCH --output=logs/SNP_STR_PVE-%j.out
#SBATCH --error=logs/SNP_STR_PVE-%j.err
#SBATCH --time=03:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=50G
#SBATCH --nodes=1
############################################

CONDITION=${1}
DATA=${2}
PC=${3}

module load R/3.6.3-no-openblas
Rscript SNP_STR_PVE_updated.R $CONDITION $DATA $PC

wait

echo " COMPLETED AT `date`"

## EOF ##
