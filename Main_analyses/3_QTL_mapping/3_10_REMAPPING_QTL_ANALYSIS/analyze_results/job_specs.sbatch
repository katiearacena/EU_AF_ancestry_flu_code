#!/bin/bash

############################################
#SBATCH --job-name=remapping_comparison
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=3:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=48G
#SBATCH --nodes=1
############################################

CONDITION=${1}
DATA=${2}
PC=${3}

module load R
Rscript compare_remapped_with_orig_results_plotting_and_corr.V2.R $CONDITION $DATA $PC

wait

echo " COMPLETED AT `date`"

## EOF ##
