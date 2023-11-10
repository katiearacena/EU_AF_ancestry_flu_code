#!/bin/bash

############################################
#SBATCH --job-name=MATRIXeQTL
#SBATCH --output=logs/ME-%j.out
#SBATCH --error=logs/ME-%j.err
#SBATCH --time=36:00:00
#SBATCH --partition=bigmem2
#SBATCH --mem=80G
#SBATCH --nodes=1
############################################

CONDITION=${1}
PC=${2}
TEMP_DIR=${3}
OUT_DIR=${4}

module load R
Rscript eQTL_mapping_MAIN.R $CONDITION $PC $TEMP_DIR $OUT_DIR

wait

echo " COMPLETED AT `date`"

## EOF ##
