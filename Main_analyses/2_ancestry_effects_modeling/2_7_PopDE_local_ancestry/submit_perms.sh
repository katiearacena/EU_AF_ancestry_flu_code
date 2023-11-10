#!/bin/bash

DATA_TYPES=(ATACseq ATACseq ATACseq ATACseq ATACseq)
PERMS=(1 2 3 4 5)


# make logs directory.
mkdir logs

LEN=${#DATA_TYPES[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       DATA=${DATA_TYPES[$NUM]}
       PERM=${PERMS[$NUM]}

       echo " ************************************* "
       echo " BEGIN POPDE PIPELINE FOR: $DATA $PERM "
       echo " ************************************* "
       sbatch perm_job_specs.sbatch $DATA $PERM&
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
