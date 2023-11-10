#!/bin/bash

############################################
#SBATCH --job-name=predictDB_model
#SBATCH --output=logs/model.%j.out
#SBATCH --error=logs/model.%j.err
#SBATCH --time=16:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=8G
#SBATCH --nodes=1
############################################



study=${1}
chrom=${2}
alpha=${3}
window=${4}
condition=${5}

mkdir logs

module load R/3.6.3-no-openblas
echo $study $chrom $condition
Rscript ./create_model.R $study $chrom $alpha $window $condition
