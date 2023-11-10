#!/bin/bash

############################################
#SBATCH --job-name=S-PrediXcan
#SBATCH --output=logs/S-PrediXcan-%j_WGBS.out
#SBATCH --error=logs/S-PrediXcan-%j_WGBS.err
#SBATCH --time=96:00:00
#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --mem=490G
#SBATCH --nodes=1
############################################




DATA=${1}
CONDITION=${2}


module load python

###### Run S-PrediXcan ######
GWAS_FILE="ADULT1_ADULT2_ONSET_ASTHMA.20180716.allchr.assoc.GC"

echo $DATA $CONDITION $GWAS_FILE 
../software/SPrediXcan.py \
--model_db_path ../../PredictDBPipeline/data/output/dbs/${DATA}_${CONDITION}_ME_snps_alpha0.5_window1e6_filtered.db \
--covariance ../../PredictDBPipeline/data/output/allCovariances/${DATA}_${CONDITION}_ME_snps_alpha0.5_window1e6.txt.gz \
--keep_non_rsid \
--gwas_folder ../GWAS \
--gwas_file_pattern ${GWAS_FILE} \
--snp_column SNP \
--effect_allele_column ALLELE1 \
--non_effect_allele_column ALLELE0 \
--beta_column BETA \
--pvalue_column P_BOLT_LMM_INF \
--output_file ../results/${GWAS_FILE}/${DATA}_${CONDITION}.csv


echo " COMPLETED AT `date`"

## EOF ##
