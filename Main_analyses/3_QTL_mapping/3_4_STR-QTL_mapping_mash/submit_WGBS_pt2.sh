#!/bin/bash
mkdir logs

PARTS=(1 2 3 4 5 6)

LEN=${#PARTS[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       PART=${PARTS[$NUM]}

       echo " ************************************** "
       echo " BEGIN WGBS MASH PIPELINE FOR: $PART "
       echo " ************************************** "
       sbatch WGBS_job_specs.pt2.sbatch $PART &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
