#!/bin/bash

#SBATCH --job-name=ldsc_continuous_annot_1KG_AFR
#SBATCH --output=ldsc_continuous_annot_1KG_AFR_%J.out
#SBATCH --error=ldsc_continuous_annot_1KG_AFR_%J.err
#SBATCH --partition=broadwl
#SBATCH --account=pi-xinhe
#SBATCH --mem=20G

## Compute annotation-specific LD scores for African population
## https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial

ANNOT_PREFIX=$1
ANNOT_DIR=$2
OUT_DIR=$3

echo "ANNOT_PREFIX: ${ANNOT_PREFIX}"
echo "ANNOT_DIR: ${ANNOT_DIR}"
echo "OUT_DIR: ${OUT_DIR}"

LDSR_DIR=/project2/xinhe/kevinluo/ldsc/

source activate ldsc

module load R/4.2.0

cd ~/projects/flu_regulatoryQTLs/code

mkdir -p ${OUT_DIR}/LDSCORE_1000G_AFR_Phase3/${ANNOT_PREFIX}_ldscores

for chrom in {1..22}
do
  echo ${chrom}

  ## Step 1: Creating ldsc-friendly annot file
  echo "Make ldsc-friendly annotation files for ${ANNOT_PREFIX}.bed"

  Rscript ../code/make_ldsc_continuous_annot.R \
  --bed ${ANNOT_DIR}/annot_bed/${ANNOT_PREFIX}.bed \
  --bim ${LDSR_DIR}/LDSCORE/1000G_AFR_Phase3_plink/1000G.AFR.phase3.${chrom}.bim \
  --out ${OUT_DIR}/LDSCORE_1000G_AFR_Phase3/${ANNOT_PREFIX}_ldscores/${ANNOT_PREFIX}.${chrom}.annot.gz \
  --out_format "full-annot"

  ## Step 2: Computing LD scores using the annot file
  echo "Computing LD scores with the annot file ${ANNOT_PREFIX}.${chrom}.annot.gz"

  python /project2/xinhe/kevinluo/software/ldsc/ldsc.py \
  --l2 \
  --bfile ${LDSR_DIR}/LDSCORE/1000G_AFR_Phase3_plink/1000G.AFR.phase3.${chrom} \
  --print-snps ${LDSR_DIR}/LDSCORE/listHM3.txt \
  --ld-wind-cm 1 \
  --annot ${OUT_DIR}/LDSCORE_1000G_AFR_Phase3/${ANNOT_PREFIX}_ldscores/${ANNOT_PREFIX}.${chrom}.annot.gz \
  --out ${OUT_DIR}/LDSCORE_1000G_AFR_Phase3/${ANNOT_PREFIX}_ldscores/${ANNOT_PREFIX}.${chrom}

done
