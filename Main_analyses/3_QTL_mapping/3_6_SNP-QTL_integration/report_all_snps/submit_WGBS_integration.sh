#!/bin/bash
mkdir logs

DATA=(7)
CONDITIONS=(NI Flu)

LEN=${#DATA[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       DAT=${DATA[$NUM]}
       CON=${CONDITIONS[$NUM]}

       echo " ************************************** "
       echo " BEGIN QTL INTEGRATION FOR: $DAT $CON "
       echo " ************************************** "
       sbatch WGBS_integration_job_specs.sbatch $DAT $CON &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
