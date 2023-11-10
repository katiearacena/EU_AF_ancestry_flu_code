#!/bin/bash

#SBATCH --job-name=torus
#SBATCH --output=torus_%J.out
#SBATCH --account=pi-xinhe
#SBATCH --partition=broadwl
#SBATCH --mem=30G

## parameters
EXP_NAME=$1

echo "Run Torus to get SNP-level prior for: ${EXP_NAME}"

module load R/4.0.4

#Rscript ~/projects/flu_regulatoryQTLs/code/get_torus_dist_enrichment.R \
#  --phenotype ${EXP_NAME} \
#  --outdir /scratch/midway2/kaixuan/flu_regulatoryQTLs/torus

Rscript ~/projects/flu_regulatoryQTLs/code/get_torus_dist_priors.R \
  --phenotype ${EXP_NAME} \
  --outdir /scratch/midway2/kaixuan/flu_regulatoryQTLs/torus

echo "Done Torus."

