#!/bin/bash

###########################################
#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=490G
#SBATCH --job-name="S-PrediXcan"
#SBATCH --output=logs/S-PrediXcan-%j_WGBS.out
#SBATCH --error=logs/S-PrediXcan-%j_WGBS.err
###########################################




DATA=${1}
CONDITION=${2}


module load python

###### Run S-PrediXcan ######
GWAS_FILE="sle_Vyse_chr1-22.txt"

echo $DATA $CONDITION $GWAS_FILE
../software/SPrediXcan.py \
--model_db_path ../../PredictDBPipeline/data/output/dbs/${DATA}_${CONDITION}_ME_snps_alpha0.5_window1e6_filtered.db \
--covariance ../../PredictDBPipeline/data/output/allCovariances/${DATA}_${CONDITION}_ME_snps_alpha0.5_window1e6.txt.gz \
--keep_non_rsid \
--gwas_folder ../GWAS \
--gwas_file_pattern ${GWAS_FILE} \
--snp_column SNP \
--effect_allele_column AlleleA \
--non_effect_allele_column AlleleB \
--beta_column Beta1Add \
--pvalue_column PValueAdd \
--output_file ../results/${GWAS_FILE}/${DATA}_${CONDITION}.csv


echo " COMPLETED AT `date`"

## EOF ##
