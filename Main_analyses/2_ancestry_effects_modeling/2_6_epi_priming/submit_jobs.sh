#!/bin/bash

DATA_TYPES=(ATACseq H3K27ac H3K27me3 H3K4me1 H3K4me3)

# make logs directory.
mkdir logs

LEN=${#DATA_TYPES[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       DATA=${DATA_TYPES[$NUM]}

       echo " ******************************* "
       echo " BEGIN ENR PIPELINE FOR: $DATA "
       echo " ******************************* "
       sbatch job_specs.sbatch $DATA &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
