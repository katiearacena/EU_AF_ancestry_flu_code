#!/bin/bash

############################################
#SBATCH --job-name=get_pred_effects
#SBATCH --output=logs/pred_effects-%j.out
#SBATCH --error=logs/pred_effects-%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --nodes=1
############################################

## WGBS needs 100 GB , all other data types 48G memory
DATA=${1}

module load R
Rscript get_predicted_pop_diffs_SNPs_only.R $DATA 

wait

echo " COMPLETED AT `date`"

## EOF ##
