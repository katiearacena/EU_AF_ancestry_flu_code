#!/bin/bash

############################################
#SBATCH --job-name=WGBS_cisregPopDE
#SBATCH --output=logs/cisreg_admixturePVE_WGBS_process_input-%j.out
#SBATCH --error=logs/cisreg_admixturePVE_WGBS_process_input-%j.err
#SBATCH --time=96:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=200G
#SBATCH --nodes=1
############################################

# make logs directory.
mkdir logs

module load R/3.6.3-no-openblas
Rscript cisreg_admixturePVE_WGBS_process_input.R

wait

echo " COMPLETED AT `date`"
