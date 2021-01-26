#!/bin/bash

############################################
#SBATCH --job-name=relaimpo
#SBATCH --output=logs/relaimpo-%j.out
#SBATCH --error=logs/relaimpo-%j.err
#SBATCH --time=96:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=128G
#SBATCH --nodes=1
#SBATCH --mail-user=karacena@uchicago.edu
#SBATCH --mail-type=ALL
############################################

Rscript pt2.R
#Rscript PVE.R
wait

echo " COMPLETED AT `date`"
