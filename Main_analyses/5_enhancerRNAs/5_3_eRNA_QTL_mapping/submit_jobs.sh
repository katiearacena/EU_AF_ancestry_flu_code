#!/bin/bash
mkdir logs


CONDITIONS=(NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu)
PCS=(0 0 1 1 2 2 3 3 4 4 5 5 6 6)
TEMP_DIRS=(0expPCs_NI 0expPCs_Flu 1expPCs_NI 1expPCs_Flu 2expPCs_NI 2expPCs_Flu 3expPCs_NI 3expPCs_Flu 4expPCs_NI 4expPCs_Flu 5expPCs_NI 5expPCs_Flu 6expPCs_NI 6expPCs_Flu)
OUT_DIRS=(0expPCs_NI 0expPCs_Flu 1expPCs_NI 1expPCs_Flu 2expPCs_NI 2expPCs_Flu 3expPCs_NI 3expPCs_Flu 4expPCs_NI 4expPCs_Flu 5expPCs_NI 5expPCs_Flu 6expPCs_NI 6expPCs_Flu)

LEN=${#CONDITIONS[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       CONDITION=${CONDITIONS[$NUM]}
       PC=${PCS[$NUM]}
       TEMP_DIR=${TEMP_DIRS[$NUM]}
       OUT_DIR=${OUT_DIRS[$NUM]}


       echo " ************************************** "
       echo " BEGIN MATRIXEQTL PIPELINE FOR: eRNA "
       echo " CONDITION = $CONDITION;"
       echo " EXP PC TO REGRESS = $PC;"
       echo " TEMP DIRECTORY = $TEMP_DIR;"
       echo " OUTPUT DIRECTORY = $OUT_DIR;"
       echo " ************************************** "
       sbatch job_specs.sbatch $CONDITION $PC $TEMP_DIR $OUT_DIR &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
