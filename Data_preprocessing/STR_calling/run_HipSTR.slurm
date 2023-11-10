#!/bin/sh
#SBATCH --partition=broadwl
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --job-name="hipSTR"
#SBATCH --output=EU_AF_hipSTR.log

echo "SLURM_JOB_ID="$SLURM_JOB_ID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR
#

#
filenames=`cat bamfilelist` 
# bamfilelist contains the names of the bam files of the samples in this project

../HipSTR --bams $filenames\
        --fasta <GRCh37_FASTA> \
        --regions GRCh37.hipstr_reference.bed.gz \
	--str-vcf str_calls_EU_AF.vcf.gz


echo "done"
#
