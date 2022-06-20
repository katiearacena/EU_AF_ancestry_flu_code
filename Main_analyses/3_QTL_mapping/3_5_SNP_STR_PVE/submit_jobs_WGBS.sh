#!/bin/bash


CONDITIONS=(NI Flu)
DATA=(WGBS WGBS)
PCS=(0 0)


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


