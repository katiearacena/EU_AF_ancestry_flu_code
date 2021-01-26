#!/bin/bash

############################################
#SBATCH --job-name=popDE_jobs
#SBATCH --output=logs/PopDE-%j.out
#SBATCH --error=logs/PopDE-%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=24G
#SBATCH --nodes=1
#SBATCH --mail-user=karacena@uchicago.edu
#SBATCH --mail-type=ALL
############################################

DATA=${1}

module load R

Rscript PopDE.R $DATA FALSE

wait

echo " COMPLETED AT `date`"

## EOF ##
