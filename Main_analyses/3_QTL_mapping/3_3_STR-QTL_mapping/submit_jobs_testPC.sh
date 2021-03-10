#!/bin/bash
mkdir logs




for PC in $(seq 0 0)
    do
       DAT=RNAseq
       CONDITION=NI
       TEMP_DIR=RNAseq_NI
       OUT_DIR=RNAseq_NI_$PC


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


for PC in $(seq 0 0)
    do
       DAT=RNAseq
       CONDITION=Flu        
       TEMP_DIR=RNAseq_Flu
       OUT_DIR=RNAseq_Flu_$PC


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




