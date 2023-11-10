#!/bin/bash

############################################
#SBATCH --job-name=eRNA_logreg_pt1
#SBATCH --output=logs/eRNA_logreg_pt1-%j.out
#SBATCH --error=logs/eRNA_logreg_pt1-%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --account=pi-\lbarreiro
############################################

mkdir logs 
module load R

Rscript QTL_eRNA_top_snps.R NI

wait

echo " COMPLETED AT `date`"

## EOF ##
