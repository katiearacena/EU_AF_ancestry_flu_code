#!/bin/bash

############################################
#SBATCH --job-name=cisregPopDE
#SBATCH --output=logs/cisSNP_STRregPopDE-%j.out
#SBATCH --error=logs/cisSNP_STRregPopDE-%j.err
#SBATCH --time=168:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=60G
#SBATCH --nodes=1
############################################

DATA=${1}

module load R/3.6.3-no-openblas

Rscript C-popDE_cisreg_SNPsSTRs.R $DATA

wait

echo " COMPLETED AT `date`"

## EOF ##
