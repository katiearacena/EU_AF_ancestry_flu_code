#!/bin/bash

############################################
#SBATCH --job-name=mash_QTL
#SBATCH --output=logs/mash_QTL-%j.out
#SBATCH --error=logs/mash_QTL-%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem2
#SBATCH --mem=300G
#SBATCH --nodes=1
#SBATCH --mail-user=karacena@uchicago.edu
#SBATCH --mail-type=ALL
############################################

DATA=${1}

module load R

Rscript mash_matrixeqtl.R $DATA

wait

echo " COMPLETED AT `date`"

## EOF ##
