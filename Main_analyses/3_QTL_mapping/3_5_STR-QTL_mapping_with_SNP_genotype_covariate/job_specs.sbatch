#!/bin/bash

############################################
#SBATCH --job-name=MATRIXeQTL
#SBATCH --output=logs/ME-%j.out
#SBATCH --error=logs/ME-%j.err
#SBATCH --time=36:00:00
#SBATCH --partition=lbarreiro
#SBATCH --mem=60G
#SBATCH --nodes=1
############################################

CONDITION=${1}
DATA=${2}
PC=${3}
TEMP_DIR=${4}
OUT_DIR=${5}

module load R/3.6.3-no-openblas
Rscript eQTL_mapping_MASTER.R $CONDITION $DATA $PC $TEMP_DIR $OUT_DIR

wait

echo " COMPLETED AT `date`"

## EOF ##
