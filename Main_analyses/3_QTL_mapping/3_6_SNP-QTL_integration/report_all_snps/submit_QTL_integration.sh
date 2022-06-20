#!/bin/bash
mkdir logs

DATA=(1 2 3 4 5 6 1 2 3 4 5 6)
CONDITIONS=(NI NI NI NI NI NI Flu Flu Flu Flu Flu Flu)

LEN=${#DATA[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       DAT=${DATA[$NUM]}
       CON=${CONDITIONS[$NUM]}

       echo " ************************************** "
       echo " BEGIN QTL INTEGRATION FOR: $DAT $CON "
       echo " ************************************** "
       sbatch integration_job_specs.sbatch $DAT $CON &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
