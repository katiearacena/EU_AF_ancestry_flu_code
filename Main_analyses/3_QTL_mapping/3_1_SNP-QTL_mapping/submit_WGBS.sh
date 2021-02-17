
#!/bin/bash
mkdir logs

CONDITIONS=(NI Flu NI Flu)
PCS=(0 0 1 1)
TEMP_DIRS=(qqnorm_0expPC_NI qqnorm_0expPC_Flu qqnorm_1expPC_NI qqnorm_1expPC_Flu)
OUT_DIRS=(qqnorm_0expPCs qqnorm_0expPCs qqnorm_1expPCs qqnorm_1expPCs)

#CONDITIONS=(NI Flu NI Flu)
#PCS=(0 0 1 1)
#TEMP_DIRS=(0expPC_NI 0expPC_Flu 1expPC_NI 1expPC_Flu)
#OUT_DIRS=(0expPCs 0expPCs 1expPCs 1exPCs)

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
