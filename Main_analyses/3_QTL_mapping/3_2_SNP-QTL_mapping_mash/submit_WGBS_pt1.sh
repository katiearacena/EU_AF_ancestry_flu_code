#!/bin/bash

## part 1 of the WGBS mash pipeline extracts and saves the data in the proper format to run mash.
## part 2 of the pipeline applies mash to the data and saves results.

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
