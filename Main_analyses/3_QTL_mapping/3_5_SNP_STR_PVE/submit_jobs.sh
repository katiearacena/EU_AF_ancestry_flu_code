#!/bin/bash

CONDITIONS=(Flu NI Flu NI Flu NI Flu NI Flu NI Flu)
DATA=(RNAseq ATACseq ATACseq H3K27ac H3K27ac H3K27me3 H3K27me3 H3K4me1 H3K4me1 H3K4me3 H3K4me3)
PCS=(7 5 5 3 2 3 3 4 5 3 4)

#CONDITIONS=(NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu)
#DATA=(RNAseq RNAseq ATACseq ATACseq H3K27ac H3K27ac H3K27me3 H3K27me3 H3K4me1 H3K4me1 H3K4me3 H3K4me3 WGBS WGBS)
#PCS=(3 7 5 5 3 2 3 3 4 5 3 4 0 0)


LEN=${#CONDITIONS[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       DAT=${DATA[$NUM]}
       CONDITION=${CONDITIONS[$NUM]}
       PC=${PCS[$NUM]}

       echo " ************************************** "
       echo " BEGIN RELAIMPO ESTIMATION FOR: $DAT "
       echo " CONDITION = $CONDITION;"
       echo " EXP PC TO REGRESS = $PC;"
       echo " ************************************** "
       sbatch job_specs.sbatch $CONDITION $DAT $PC &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##


