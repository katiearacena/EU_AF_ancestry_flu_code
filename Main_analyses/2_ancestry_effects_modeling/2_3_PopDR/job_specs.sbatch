#!/bin/bash

############################################
#SBATCH --job-name=popDR_jobs
#SBATCH --output=logs/PopDR-%j.out
#SBATCH --error=logs/PopDR-%j.err
#SBATCH --time=10:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=24G
#SBATCH --nodes=1
############################################

DATA=${1}

module load R

Rscript PopDR.R $DATA FALSE

wait

echo " COMPLETED AT `date`"

## EOF ##
