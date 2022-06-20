#!/bin/bash

############################################
#SBATCH --job-name=mash_QTL.WGBS
#SBATCH --output=logs/mash_QTL.WGBS-%j.out
#SBATCH --error=logs/mash_QTL.WGBS-%j.err
#SBATCH --time=128:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --mail-user=ontalin@uchicago.edu
#SBATCH --mail-type=ALL
############################################


module load R/3.6.3-no-openblas

Rscript WGBS_mash_pt1.R

wait

echo " COMPLETED AT `date`"

## EOF ##
