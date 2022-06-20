#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10GB

module load python

rgt-hint footprinting --atac-seq \
--output-location=/project2/lbarreiro/users/katie/EU_AF_ancestry_flu/HINT/outputs/meta/footprints/ \
--output-prefix=$1 \
/project2/lbarreiro/users/katie/EU_AF_ancestry_flu/HINT/outputs/meta/Flu_meta.bam \
$2
