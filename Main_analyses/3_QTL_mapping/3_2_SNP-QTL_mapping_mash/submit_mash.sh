
#!/bin/bash
mkdir logs

DATA=(RNAseq ATACseq H3K27ac H3K27me3 H3K4me1 H3K4me3)

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
