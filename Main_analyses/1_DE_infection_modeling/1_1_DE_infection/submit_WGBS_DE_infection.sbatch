#!/bin/bash

############################################
#SBATCH --job-name=WGBS_InfectionDE
#SBATCH --output=logs/WGBS_InfectionDE-%j.out
#SBATCH --error=logs/WGBS_InfectionDE-%j.err
#SBATCH --time=96:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=400G
#SBATCH --nodes=1
#SBATCH --mail-user=karacena@uchicago.edu
#SBATCH --mail-type=ALL
############################################

# make logs directory.
mkdir logs 

Rscript WGBS_DE_infection.R
wait

echo " COMPLETED AT `date`"

## EOF ##
