#!/bin/bash

############################################
#SBATCH --job-name=mash_QTL.WGBS
#SBATCH --output=logs/mash_QTL.WGBS-%j.out
#SBATCH --error=logs/mash_QTL.WGBS-%j.err
#SBATCH --time=128:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=200G
#SBATCH --nodes=1
############################################

PART=${1}

module load R

Rscript WGBS_mash_pt2.R $PART

wait

echo " COMPLETED AT `date`"

## EOF ##
