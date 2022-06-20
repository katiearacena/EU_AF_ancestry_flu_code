#!/bin/bash
mkdir logs

DATA=(RNAseq ATACseq H3K27ac H3K27me3 H3K4me1 H3K4me3)

LEN=${#DATA[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       DAT=${DATA[$NUM]}

       echo " ************************************** "
       echo " BEGIN POPDE CISREG ANALYSIS FOR: $DAT "
       echo " ************************************** "

        sbatch C-job_specs_SNPsSTRs.sbatch  $DAT &
	
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
