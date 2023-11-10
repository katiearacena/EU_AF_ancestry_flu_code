#!/bin/bash

#SBATCH --job-name=sldsc
#SBATCH --output=sldsc_%J.out
#SBATCH --error=sldsc_%J.err
#SBATCH --account=pi-xinhe
#SBATCH --partition=broadwl
#SBATCH --mem=10G

GWAS_DIR=$1
TRAIT=$2
ANNOT_PREFIX=$3
ANNOT_BASELINE=$4
ANNOT_LDSCORE_DIR=$5
SLDSC_OUTDIR=$6

LDSC_DIR=/project2/xinhe/kevinluo/ldsc

BASELINELD_DIR=/project2/xinhe/kevinluo/ldsc/annot/ldscores/${ANNOT_BASELINE}

source activate ldsc

echo "GWAS_DIR: ${GWAS_DIR}"
echo "TRAIT: ${TRAIT}"
echo "ANNOT_PREFIX: ${ANNOT_PREFIX}"
echo "ANNOT_BASELINE: ${ANNOT_BASELINE}"
echo "ANNOT_LDSCORE_DIR: ${ANNOT_LDSCORE_DIR}"
echo "SLDSC_OUTDIR: ${SLDSC_OUTDIR}"
echo "BASELINELD_DIR: ${BASELINELD_DIR}"

out_dir=${SLDSC_OUTDIR}/${TRAIT}/${ANNOT_BASELINE}
mkdir -p ${out_dir}

if [ ! -d "${ANNOT_LDSCORE_DIR}" ]; then
  echo "Error: Directory ${ANNOT_LDSCORE_DIR} not found."
  exit 1
fi

python /project2/xinhe/kevinluo/software/ldsc/ldsc.py \
--h2 ${GWAS_DIR}/${TRAIT}.sumstats.gz \
--ref-ld-chr ${BASELINELD_DIR}/${ANNOT_BASELINE}.,${ANNOT_LDSCORE_DIR}/${ANNOT_PREFIX}_ldscores/${ANNOT_PREFIX}. \
--frqfile-chr ${LDSC_DIR}/LDSCORE/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ${LDSC_DIR}/LDSCORE/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot --print-cov --print-coefficients --print-delete-vals \
--out ${out_dir}/${TRAIT}_${ANNOT_PREFIX}_${ANNOT_BASELINE}
