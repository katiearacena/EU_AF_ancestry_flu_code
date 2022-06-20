
#!/bin/bash
mkdir logs

CONDITIONS=(NI Flu)
PCS=(0 0)
TEMP_DIRS=(WGBS_NI WGBS_Flu)
OUT_DIRS=(WGBS_NI WGBS_Flu)

LEN=${#CONDITIONS[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do

       CONDITION=${CONDITIONS[$NUM]}
       PC=${PCS[$NUM]}
       TEMP_DIR=${TEMP_DIRS[$NUM]}
       OUT_DIR=${OUT_DIRS[$NUM]}


       echo " ************************************** "
       echo " BEGIN MATRIXEQTL PIPELINE FOR: WGBS "
       echo " CONDITION = $CONDITION;"
       echo " EXP PC TO REGRESS = $PC;"
       echo " TEMP DIRECTORY = $TEMP_DIR;"
       echo " OUTPUT DIRECTORY = $OUT_DIR;"
       echo " ************************************** "
       sbatch WGBS_job_specs.sbatch $CONDITION $PC $TEMP_DIR $OUT_DIR &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
