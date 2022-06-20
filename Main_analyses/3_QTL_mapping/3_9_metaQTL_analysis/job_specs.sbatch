#!/bin/bash

############################################
#SBATCH --job-name=EP
#SBATCH --output=logs/EP-%j.out
#SBATCH --error=logs/EP-%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=48G
#SBATCH --nodes=1
############################################

DATA=${1}

Rscript epi_NI_levels_associated_with_caQTL.R $DATA

wait

echo " COMPLETED AT `date`"

## EOF ##
