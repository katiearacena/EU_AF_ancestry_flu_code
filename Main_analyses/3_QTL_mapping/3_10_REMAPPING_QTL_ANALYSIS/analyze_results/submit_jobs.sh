#!/bin/bash
mkdir logs

#CONDITIONS=(NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu)
#DATA=(RNAseq RNAseq ATACseq ATACseq H3K27ac H3K27ac H3K27me3 H3K27me3 H3K4me1 H3K4me1 H3K4me3 H3K4me3)
#PCS=(4 4 3 3 1 1 1 2 2 4 2 2)

CONDITIONS=(NI Flu)
DATA=(ATACseq ATACseq)
PCS=(3 3)

LEN=${#CONDITIONS[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       DAT=${DATA[$NUM]}
       CONDITION=${CONDITIONS[$NUM]}
       PC=${PCS[$NUM]}
       echo " ************************************** "
       echo " BEGIN COMPARISON FOR: $DAT "
       echo " CONDITION = $CONDITION;"
       echo " EXP PC = $PC;"
       echo " ************************************** "
       sbatch job_specs.sbatch $CONDITION $DAT $PC &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
