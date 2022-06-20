#!/bin/bash
mkdir logs

DATA=(RNAseq ATACseq H3K27ac H3K27me3 H3K4me1 H3K4me3 WGBS RNAseq ATACseq H3K27ac H3K27me3 H3K4me1 H3K4me3 WGBS)
CONDITIONS=(NI NI NI NI NI NI NI Flu Flu Flu Flu Flu Flu Flu)

LEN=${#DATA[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       DAT=${DATA[$NUM]}
       CON=${CONDITIONS[$NUM]}

       echo " ************************************** "
       echo " BEGIN TF ENRICHMENT PIPELINE FOR: $DAT $CON "
       echo " ************************************** "
       sbatch job_specs.sbatch $DAT $CON &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
