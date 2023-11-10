#!/bin/bash
mkdir logs


DATAs="RNAseq ATACseq H3K27ac H3K27me3 H3K4me1 H3K4me3"


CONDITION="Flu"
for DATA in $DATAs
    do
       echo " ************************************** "
       echo " BEGIN PREDIXCAN FOR: $DATA $CONDITION"
       echo " ************************************** "
       sbatch job_specs.sbatch $DATA $CONDITION
    done
wait


CONDITION="NI"  
for DATA in $DATAs
    do
       echo " ************************************** "
       echo " BEGIN PREDIXCAN FOR: $DATA $CONDITION"
       echo " ************************************** "
       sbatch job_specs.sbatch $DATA $CONDITION
    done
wait



echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##


