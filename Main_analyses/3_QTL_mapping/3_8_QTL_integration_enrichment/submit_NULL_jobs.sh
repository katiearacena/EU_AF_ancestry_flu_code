#!/bin/bash
mkdir logs

## run for each chromosome and condition separtely
DATA=(1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 16 17 17 18 18 19 19 20 20 21 21 22 22)
CONDITIONS=(NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu NI Flu)


LEN=${#DATA[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       DAT=${DATA[$NUM]}
       CON=${CONDITIONS[$NUM]}

       echo " ************************************** "
       echo " BEGIN TF ENRICHMENT PIPELINE FOR: $DAT $CON "
       echo " ************************************** "
       sbatch QTL_integ_footprinting_NULL_job_specs.sbatch $DAT $CON &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
