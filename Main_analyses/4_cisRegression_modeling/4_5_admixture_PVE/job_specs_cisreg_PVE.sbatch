#!/bin/bash

############################################
#SBATCH --job-name=cisregPopDE
#SBATCH --output=logs/cisreg_PVE-%j.out
#SBATCH --error=logs/cisreg_PVE-%j.err
#SBATCH --time=168:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=60G
#SBATCH --nodes=1
############################################

DATA=${1}

module load R/3.6.3-no-openblas

Rscript cisreg_admixturePVE.R $DATA

wait

echo " COMPLETED AT `date`"

## EOF ##
