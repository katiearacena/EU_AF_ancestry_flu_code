#!/bin/bash

############################################
#SBATCH --job-name=SAMEpred_effects
#SBATCH --output=logs/SAMEeffects-%j.out
#SBATCH --error=logs/SAMEeffects-%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=100G
#SBATCH --nodes=1
############################################

CONDITION=${1}
DATA=${2}
PC=${3}

module load R/3.6.3-no-openblas
Rscript get_predicted_pop_diffs_SNPs_STRs_same_model.R $CONDITION $DATA $PC

wait

echo " COMPLETED AT `date`"

## EOF ##
