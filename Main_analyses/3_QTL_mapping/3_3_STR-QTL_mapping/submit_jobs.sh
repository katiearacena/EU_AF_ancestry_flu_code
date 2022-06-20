
#!/bin/bash
mkdir logs


CONDITIONS=(NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu)
DATA=(RNAseq RNAseq ATACseq ATACseq H3K27ac H3K27ac H3K27me3 H3K27me3 H3K4me1 H3K4me1 H3K4me3 H3K4me3)
PCS=(3 7 5 5 3 2 3 3 4 5 3 4)
TEMP_DIRS=(RNAseq_NI RNAseq_Flu ATACseq_NI ATACseq_Flu H3K27ac_NI H3K27ac_Flu H3K27me3_NI H3K27me3_Flu H3K4me1_NI H3K4me1_Flu H3K4me3_NI H3K4me3_Flu)
OUT_DIRS=(RNAseq_NI RNAseq_Flu ATACseq_NI ATACseq_Flu H3K27ac_NI H3K27ac_Flu H3K27me3_NI H3K27me3_Flu H3K4me1_NI H3K4me1_Flu H3K4me3_NI H3K4me3_Flu)

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

