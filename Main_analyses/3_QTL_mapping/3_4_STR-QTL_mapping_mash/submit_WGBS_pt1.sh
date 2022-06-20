#!/bin/bash
mkdir logs

LEN=${#PARTS[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       PART=${PARTS[$NUM]}

       echo " ************************************** "
       echo " BEGIN MASH PIPELINE FOR: WGBS pt 1 "
       echo " ************************************** "
       sbatch WGBS_job_specs.pt1.sbatch &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
