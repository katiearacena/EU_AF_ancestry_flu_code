
#!/bin/bash
mkdir logs

DATA=(H3K27ac)

LEN=${#DATA[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       DAT=${DATA[$NUM]}

       echo " ************************************** "
       echo " BEGIN MASH PIPELINE FOR: $DAT "
       echo " ************************************** "
       sbatch job_specs.sbatch $DAT &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
