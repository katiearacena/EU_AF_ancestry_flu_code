#!/bin/bash

############################################
#SBATCH --job-name=GSEA
#SBATCH --output=logs/GSEA-%j.out
#SBATCH --error=logs/GSEA-%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=48G
#SBATCH --nodes=1
#SBATCH --mail-user=karacena@uchicago.edu
#SBATCH --mail-type=ALL
############################################

Rscript infection_GSEA.R
wait

echo " COMPLETED AT `date`"
