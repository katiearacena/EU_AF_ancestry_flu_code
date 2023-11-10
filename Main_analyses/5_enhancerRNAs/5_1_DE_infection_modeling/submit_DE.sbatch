#!/bin/bash

# This is used by submit_DE_infection.sh to indicate the partition information.

############################################
#SBATCH --job-name=InfectionDE
#SBATCH --output=logs/InfectionDE-%j.out
#SBATCH --error=logs/InfectionDE-%j.err
#SBATCH --time=4:00:00
#SBATCH --partition=lbarreiro
#SBATCH --account=pi-\lbarreiro
#SBATCH --mem=16G
#SBATCH --nodes=1
############################################

mkdir logs

module load R

Rscript filter_DE_infection_effects.R

wait

echo " COMPLETED AT `date`"

## EOF ##
