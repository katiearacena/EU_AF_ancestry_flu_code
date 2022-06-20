#!/bin/bash

############################################
#SBATCH --job-name=QTL_perm
#SBATCH --output=logs/QTL_perm-%j.out
#SBATCH --error=logs/QTL_perm-%j.err
#SBATCH --time=168:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=200G
#SBATCH --nodes=1
############################################

mkdir logs 

DATA=${1}
CONDITION=${2}

Rscript QTL_permutation.R $DATA $CONDITION

wait

echo " COMPLETED AT `date`"

## EOF ##
