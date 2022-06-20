#!/bin/bash

############################################
#SBATCH --job-name=WGBS_MATRIXeQTL
#SBATCH --output=logs/WGBS_ME-%j.out
#SBATCH --error=logs/WGBS_ME-%j.err
#SBATCH --time=96:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=300G
#SBATCH --nodes=1
############################################

CONDITION=${1}
PC=${2}
TEMP_DIR=${3}
OUT_DIR=${4}

module load R

Rscript WGBS_qqnorm_matrixeqtl.R $CONDITION $PC $TEMP_DIR $OUT_DIR

wait

echo " COMPLETED AT `date`"

## EOF ##
