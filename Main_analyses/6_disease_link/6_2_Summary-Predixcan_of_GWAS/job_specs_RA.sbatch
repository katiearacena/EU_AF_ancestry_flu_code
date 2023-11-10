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
GWAS_FILE="RA_GWASmeta_TransEthnic_v2.txt"

echo $DATA $CONDITION $GWAS_FILE
../software/SPrediXcan.py \
--model_db_path ../../PredictDBPipeline/data/output/dbs/${DATA}_${CONDITION}_ME_snps_alpha0.5_window1e6_filtered.db \
--covariance ../../PredictDBPipeline/data/output/allCovariances/${DATA}_${CONDITION}_ME_snps_alpha0.5_window1e6.txt.gz \
--keep_non_rsid \
--gwas_folder ../GWAS \
--gwas_file_pattern ${GWAS_FILE} \
--snp_column SNP \
--effect_allele_column A1 \
--non_effect_allele_column A2 \
--or_column 'OR(A1)' \
--pvalue_column P-val \
--output_file ../results/${GWAS_FILE}/${DATA}_${CONDITION}.csv


echo " COMPLETED AT `date`"

## EOF ##
