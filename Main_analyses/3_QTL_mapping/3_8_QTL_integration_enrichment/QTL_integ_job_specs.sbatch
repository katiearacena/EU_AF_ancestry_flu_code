#!/bin/bash

############################################
#SBATCH --job-name=TF_enrich
#SBATCH --output=logs/TF_enrich-%j.out
#SBATCH --error=logs/TF_enrich-%j.err
#SBATCH --time=128:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=48G
#SBATCH --nodes=1
############################################

QTLNUM=${1}
CONDITION=${2}

module load R

Rscript QTL_integ_TF_footprinting_analysis.R $QTLNUM $CONDITION

wait
echo " COMPLETED AT `date`"
## EOF ##
