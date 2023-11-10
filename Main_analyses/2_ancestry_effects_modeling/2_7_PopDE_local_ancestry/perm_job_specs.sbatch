#!/bin/bash

############################################
#SBATCH --job-name=LA_popDE_perm
#SBATCH --output=logs/LA_PopDE-%j.out
#SBATCH --error=logs/LA_PopDE-%j.err
#SBATCH --account=pi-\lbarreiro
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --mem=24G
#SBATCH --nodes=1
############################################

DATA=${1}
PERM_NUM=${2}

module load gsl
module load R

Rscript PopDE_LA.perm.R $DATA $PERM_NUM

wait

echo " COMPLETED AT `date`"

## EOF ##
