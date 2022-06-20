#!/bin/bash

# make logs directory.
mkdir logs

MODELS=(popDR_RNAseq_IFNA_score popDR_RNAseq_IFNG_score popDR_RNAseq_TNFA_score popDR_RNAseq_IL6_score)
ALPHAS=(0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1)

LEN1=${#MODELS[@]}
LEN2=${#ALPHAS[@]}


for (( NUM1=0; NUM1<$LEN1; NUM1++ ))
	do
		MODEL=${MODELS[$NUM1]}
		
		for (( NUM2=0; NUM2<$LEN2; NUM2++ ))
    		do
       			ALPHA=${ALPHAS[$NUM2]}

       			echo " ********************************************* "
       			echo " BEGIN ENR PIPELINE FOR: $MODEL "alpha="$ALPHA "
       			echo " ********************************************* "
                sbatch WGBS_job_specs.sbatch $MODEL $ALPHA &
       		done
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
