#!/bin/bash

## part 1 of the WGBS mash pipeline extracts and saves the data in the proper format to run mash.
## part 2 of the pipeline applies mash to the data and saves results.

mkdir logs

PARTS=(1 2 3 4 5 6 7 8)

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
