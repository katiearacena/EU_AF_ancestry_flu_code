#!/bin/bash
mkdir logs


#CONDITIONS=(NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu)
#DATA=(RNAseq RNAseq ATACseq ATACseq H3K27ac H3K27ac H3K27me3 H3K27me3 H3K4me1 H3K4me1 H3K4me3 H3K4me3)
#PCS=(4 4 3 3 1 1 1 2 2 4 2 2)
#TEMP_DIRS=(RNAseq_NI RNAseq_Flu ATACseq_NI ATACseq_Flu H3K27ac_NI H3K27ac_Flu H3K27me3_NI H3K27me3_Flu H3K4me1_NI H3K4me1_Flu H3K4me3_NI H3K4me3_Flu)
#OUT_DIRS=(RNAseq_NI RNAseq_Flu ATACseq_NI ATACseq_Flu H3K27ac_NI H3K27ac_Flu H3K27me3_NI H3K27me3_Flu H3K4me1_NI H3K4me1_Flu H3K4me3_NI H3K4me3_Flu)

CONDITIONS=(NI Flu NI Flu NI Flu)
DATA=(RNAseq RNAseq H3K27me3 H3K27me3 H3K4me3 H3K4me3)
PCS=(4 4 1 2 2 2)
TEMP_DIRS=(RNAseq_NI RNAseq_Flu H3K27me3_NI H3K27me3_Flu H3K4me3_NI H3K4me3_Flu)
OUT_DIRS=(RNAseq_NI RNAseq_Flu H3K27me3_NI H3K27me3_Flu H3K4me3_NI H3K4me3_Flu)

#CONDITIONS=(NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu)
#DATA=(ATACseq ATACseq ATACseq ATACseq ATACseq ATACseq ATACseq ATACseq ATACseq ATACseq ATACseq ATACseq)
#PCS=(0 0 1 1 2 2 3 3 4 4 5 5)
#TEMP_DIRS=(ATACseq_NI_0 ATACseq_Flu_0 ATACseq_NI_1 ATACseq_Flu_1 ATACseq_NI_2 ATACseq_Flu_2 ATACseq_NI_3 ATACseq_Flu_3 ATACseq_NI_4 ATACseq_Flu_4 ATACseq_NI_5 ATACseq_Flu_5)
#OUT_DIRS=(ATACseq_NI_0 ATACseq_Flu_0 ATACseq_NI_1 ATACseq_Flu_1 ATACseq_NI_2 ATACseq_Flu_2 ATACseq_NI_3 ATACseq_Flu_3 ATACseq_NI_4 ATACseq_Flu_4 ATACseq_NI_5 ATACseq_Flu_5)

LEN=${#CONDITIONS[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       DAT=${DATA[$NUM]}
       CONDITION=${CONDITIONS[$NUM]}
       PC=${PCS[$NUM]}
       TEMP_DIR=${TEMP_DIRS[$NUM]}
       OUT_DIR=${OUT_DIRS[$NUM]}


       echo " ************************************** "
       echo " BEGIN MATRIXEQTL PIPELINE FOR: $DAT "
       echo " CONDITION = $CONDITION;"
       echo " EXP PC TO REGRESS = $PC;"
       echo " TEMP DIRECTORY = $TEMP_DIR;"
       echo " OUTPUT DIRECTORY = $OUT_DIR;"
       echo " ************************************** "
       sbatch job_specs.sbatch $CONDITION $DAT $PC $TEMP_DIR $OUT_DIR &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##


