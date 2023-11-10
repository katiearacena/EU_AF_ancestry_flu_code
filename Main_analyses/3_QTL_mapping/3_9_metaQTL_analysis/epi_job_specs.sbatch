#!/bin/bash

############################################
#SBATCH --job-name=epiQTL_txnresponse
#SBATCH --output=logs/epiQTL_txnresponse-%j.out
#SBATCH --error=logs/epiQTL_txnresponse-%j.err
#SBATCH --account=pi-\lbarreiro
#SBATCH --partition=caslake
#SBATCH --time=36:00:00
#SBATCH --mem=60G
#SBATCH --nodes=1
############################################

DATA=${1}


module load R
Rscript epiQTL_NI_associated_with_txn_response_connected_through_SNP.R $DATA

wait

echo " COMPLETED AT `date`"

## EOF ##
