#!/bin/bash

#SBATCH --job-name=susie_WGBS
#SBATCH --output=susie_WGBS_%J.out
#SBATCH --error=susie_WGBS_%J.err
#SBATCH --account=pi-mstephens
#SBATCH --partition=broadwl
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6

# --ntasks-per-node will be used in doParallel.R to specify the number
# of cores to use on the machine.

## parameters
EXP_NAME=$1

echo "SuSiE finemapping for: ${EXP_NAME}"

module load R/4.0.4

# Run susie with distance prior
echo "Run susie with distance prior"

Rscript ~/projects/flu_regulatoryQTLs/code/susie_allgenes_prior_WGBS.R \
  --phenotype ${EXP_NAME} \
  --L 3 \
  --susie_prior dist \
  --torus_dir /project2/xinhe/kevinluo/flu_regulatoryQTLs/torus/${EXP_NAME} \
  --susie_input_dir /scratch/midway2/kaixuan/flu_regulatoryQTLs/susie/input/${EXP_NAME} \
  --outdir /project2/xinhe/kevinluo/flu_regulatoryQTLs/susie/output/${EXP_NAME} \
  --ncores 10

# Run susie with uniform prior
echo "Run susie with uniform prior"

Rscript ~/projects/flu_regulatoryQTLs/code/susie_allgenes_prior_WGBS.R \
  --phenotype ${EXP_NAME} \
  --L 3 \
  --susie_prior uniform \
  --susie_input_dir /scratch/midway2/kaixuan/flu_regulatoryQTLs/susie/input/${EXP_NAME} \
  --outdir /project2/xinhe/kevinluo/flu_regulatoryQTLs/susie/output/${EXP_NAME} \
  --ncores 10
