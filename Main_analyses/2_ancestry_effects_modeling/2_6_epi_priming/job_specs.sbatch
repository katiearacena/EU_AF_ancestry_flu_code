#!/bin/bash

############################################
#SBATCH --job-name=epi_priming_ENR
#SBATCH --output=logs/ENR-%j.out
#SBATCH --error=logs/ENR-%j.err
#SBATCH --time=96:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=36G
#SBATCH --nodes=1
############################################

DATA=${1}

module load R

Rscript elastic_net.R $DATA

wait

echo " COMPLETED AT `date`"

## EOF ##
