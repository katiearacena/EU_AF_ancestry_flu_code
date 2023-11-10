#!/bin/bash

############################################
#SBATCH --job-name=popDE_jobs
#SBATCH --output=logs/PopDE-%j.out
#SBATCH --error=logs/PopDE-%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=lbarreiro
#SBATCH --account=pi-\lbarreiro
#SBATCH --mem=24G
#SBATCH --nodes=1
############################################

module load R
module load gsl

## run script with reproductions == FALSE
Rscript PopDE.R FALSE

wait

echo " COMPLETED AT `date`"

## EOF ##
