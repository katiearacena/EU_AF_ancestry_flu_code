#!/bin/bash

# This is used by submit_DE_infection.sh to indicate the partition information.

############################################
#SBATCH --job-name=InfectionDE
#SBATCH --output=logs/InfectionDE-%j.out
#SBATCH --error=logs/InfectionDE-%j.err
#SBATCH --time=22:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=16G
#SBATCH --nodes=1
############################################

DATA=${1}

module load R

Rscript DE_infection.R $DATA FALSE

wait

echo " COMPLETED AT `date`"

## EOF ##
