#!/bin/bash

############################################
#SBATCH --job-name=STRpred_effects
#SBATCH --output=logs/STReffects-%j.out
#SBATCH --error=logs/STReffects-%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=48G
#SBATCH --nodes=1
############################################

DATA=${1}

module load R/3.6.3-no-openblas
Rscript get_predicted_pop_diffs_STRs.R $DATA 

wait

echo " COMPLETED AT `date`"

## EOF ##
