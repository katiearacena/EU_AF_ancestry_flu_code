#!/bin/bash

############################################
#SBATCH --job-name=relaimpo
#SBATCH --output=logs/relaimpo-%j.out
#SBATCH --error=logs/relaimpo-%j.err
#SBATCH --time=96:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=128G
#SBATCH --nodes=1
############################################

module load R

Rscript PVE.R

wait

echo " COMPLETED AT `date`"
