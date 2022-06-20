#!/bin/bash
mkdir logs

DATA=(QTL1 QTL1 QTL2 QTL2 QTL3 QTL3 QTL4 QTL4 QTL5 QTL5 QTL6 QTL6 QTL7 QTL7)
CONDITIONS=(NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu)


LEN=${#DATA[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       DAT=${DATA[$NUM]}
       CON=${CONDITIONS[$NUM]}

       echo " ************************************** "
       echo " BEGIN TF ENRICHMENT PIPELINE FOR: $DAT $CON "
       echo " ************************************** "
       sbatch QTL_integ_job_specs.sbatch $DAT $CON &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
