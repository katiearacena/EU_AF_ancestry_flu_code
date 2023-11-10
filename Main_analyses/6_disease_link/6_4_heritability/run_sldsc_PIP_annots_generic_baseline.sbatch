#!/bin/bash

#SBATCH --job-name=run_sLDSC_PIP
#SBATCH --output=run_sLDSC_PIP_%J.out
#SBATCH --error=run_sLDSC_PIP_%J.err
#SBATCH --account=pi-xinhe
#SBATCH --partition=broadwl
#SBATCH --mem=20G

TRAIT=$1
ANNOT_BASELINE=$2

GWAS_DIR=/project2/xinhe/kevinluo/flu_regulatoryQTLs/GWAS_summary_stats/ldsc_format
ANNOT_LDSCORE_DIR=/project2/xinhe/kevinluo/flu_regulatoryQTLs/sldsc_results/LDSCORE
SLDSC_OUTDIR=/project2/xinhe/kevinluo/flu_regulatoryQTLs/sldsc_results/output

for bedfile in /project2/xinhe/kevinluo/flu_regulatoryQTLs/sldsc_results/annot/annot_bed/*_Flu_PIP_dist_prior_FDR0.1.bed
do
  ANNOT_PREFIX="$(basename ${bedfile} .bed)"
  echo ${ANNOT_PREFIX}
  sh ~/projects/flu_regulatoryQTLs/code/sldsc_annot_generic_baseline_separate.sbatch \
    ${GWAS_DIR} ${TRAIT} ${ANNOT_PREFIX} ${ANNOT_BASELINE} ${ANNOT_LDSCORE_DIR} ${SLDSC_OUTDIR}
done

for bedfile in /project2/xinhe/kevinluo/flu_regulatoryQTLs/sldsc_results/annot/annot_bed/*_NI_PIP_dist_prior_FDR0.1.bed
do
  ANNOT_PREFIX="$(basename ${bedfile} .bed)"
  echo ${ANNOT_PREFIX}
  sh ~/projects/flu_regulatoryQTLs/code/sldsc_annot_generic_baseline_separate.sbatch \
    ${GWAS_DIR} ${TRAIT} ${ANNOT_PREFIX} ${ANNOT_BASELINE} ${ANNOT_LDSCORE_DIR} ${SLDSC_OUTDIR}
done


#for bedfile in /project2/xinhe/kevinluo/flu_regulatoryQTLs/sldsc_results/annot/annot_bed/WGBS_Flu_full_*_PIP_dist_prior_FDR0.1.bed
#do
#  ANNOT_PREFIX="$(basename ${bedfile} .bed)"
#  echo ${ANNOT_PREFIX}
#  sh ~/projects/flu_regulatoryQTLs/code/sldsc_annot_generic_baseline_separate.sbatch \
#    ${GWAS_DIR} ${TRAIT} ${ANNOT_PREFIX} ${ANNOT_BASELINE} ${ANNOT_LDSCORE_DIR} ${SLDSC_OUTDIR}
#done

#for bedfile in /project2/xinhe/kevinluo/flu_regulatoryQTLs/sldsc_results/annot/annot_bed/WGBS_NI_full_*_PIP_dist_prior_FDR0.1.bed
#do
#  ANNOT_PREFIX="$(basename ${bedfile} .bed)"
#  echo ${ANNOT_PREFIX}
#  sh ~/projects/flu_regulatoryQTLs/code/sldsc_annot_generic_baseline_separate.sbatch \
#    ${GWAS_DIR} ${TRAIT} ${ANNOT_PREFIX} ${ANNOT_BASELINE} ${ANNOT_LDSCORE_DIR} ${SLDSC_OUTDIR}
#done

#for bedfile in /project2/xinhe/kevinluo/flu_regulatoryQTLs/sldsc_results/annot/annot_bed/*_Flu_PIP_uniform_prior_FDR0.1.bed
#do
#  ANNOT_PREFIX="$(basename ${bedfile} .bed)"
#  echo ${ANNOT_PREFIX}
#  sh ~/projects/flu_regulatoryQTLs/code/sldsc_annot_generic_baseline_separate.sbatch \
#    ${GWAS_DIR} ${TRAIT} ${ANNOT_PREFIX} ${ANNOT_BASELINE} ${ANNOT_LDSCORE_DIR} ${SLDSC_OUTDIR}
#done

#for bedfile in /project2/xinhe/kevinluo/flu_regulatoryQTLs/sldsc_results/annot/annot_bed/*_NI_PIP_uniform_prior_FDR0.1.bed
#do
# ANNOT_PREFIX="$(basename ${bedfile} .bed)"
#  echo ${ANNOT_PREFIX}
#  sh ~/projects/flu_regulatoryQTLs/code/sldsc_annot_generic_baseline_separate.sbatch \
#    ${GWAS_DIR} ${TRAIT} ${ANNOT_PREFIX} ${ANNOT_BASELINE} ${ANNOT_LDSCORE_DIR} ${SLDSC_OUTDIR}
#done
