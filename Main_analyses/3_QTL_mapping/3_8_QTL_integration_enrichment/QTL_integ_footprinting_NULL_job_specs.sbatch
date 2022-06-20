#!/bin/bash

############################################
#SBATCH --job-name=Null2_TF_enrich
#SBATCH --output=logs/Null2_TF_enrich-%j.out
#SBATCH --error=logs/Null2_TF_enrich-%j.err
#SBATCH --time=128:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=48G
#SBATCH --nodes=1
############################################

CHRNUM=${1}
CONDITION=${2}

module load R

Rscript QTL_integ_footprinting_NULL.R $CHRNUM $CONDITION

wait
echo " COMPLETED AT `date`"
## EOF ##