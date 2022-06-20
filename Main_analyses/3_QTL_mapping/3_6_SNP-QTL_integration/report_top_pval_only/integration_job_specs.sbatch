#!/bin/bash

############################################
#SBATCH --job-name=QTL_integ
#SBATCH --output=logs/QTL_integ-%j.out
#SBATCH --error=logs/QTL_integ-%j.err
#SBATCH --time=168:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=200G
#SBATCH --nodes=1
############################################

mkdir logs 

DATA=${1}
CONDITION=${2}

module load R
Rscript QTL_integration.R $DATA $CONDITION

wait

echo " COMPLETED AT `date`"

## EOF ##
