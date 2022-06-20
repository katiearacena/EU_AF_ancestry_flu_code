#!/bin/bash

############################################
#SBATCH --job-name=WGBS_epi_priming_ENR
#SBATCH --output=logs/WGBS_ENR-%j.out
#SBATCH --error=logs/WGBS_ENR-%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=100G
#SBATCH --nodes=1
############################################

MODEL_NAME=${1}
ALPHA=${2}

module load R

Rscript elastic_net_WGBS.R $MODEL_NAME $ALPHA

wait

echo " COMPLETED AT `date`"

## EOF ##
