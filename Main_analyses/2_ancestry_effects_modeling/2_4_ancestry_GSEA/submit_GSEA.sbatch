#!/bin/bash

############################################
#SBATCH --job-name=GSEA
#SBATCH --output=logs/GSEA-%j.out
#SBATCH --error=logs/GSEA-%j.err
#SBATCH --time=06:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=48G
#SBATCH --nodes=1
############################################

module load R

Rscript popDE_NI_GSEA.R
Rscript popDE_Flu_GSEA.R
Rscript popDR_GSEA.R

wait

echo " COMPLETED AT `date`"
