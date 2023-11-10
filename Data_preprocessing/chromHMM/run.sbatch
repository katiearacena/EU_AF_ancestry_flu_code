#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10GB
#SBATCH -e /project2/lbarreiro/users/katie/EU_AF_ancestry_flu/ChromHMM/logs/ChromHMM_pilot_bin%J.e
#SBATCH -o /project2/lbarreiro/users/katie/EU_AF_ancestry_flu/ChromHMM/logs/ChromHMM_pilot_bin%J.o

java -mx8000M -jar /project2/lbarreiro/programs/chromHMM/ChromHMM/ChromHMM.jar BinarizeBam -gzip /project2/lbarreiro/programs/chromHMM/ChromHMM/CHROMSIZES/hg19.txt /project2/lbarreiro/users/katie/EU_AF_ancestry_flu/ChromHMM/data/KA_CHIPseq_bams /project2/lbarreiro/users/katie/EU_AF_ancestry_flu/ChromHMM/data/KA_CHIPseq_bams/cellmarkfiletable.txt /project2/lbarreiro/users/katie/EU_AF_ancestry_flu/ChromHMM/binary_intermediates/KA_CHIPseq_bams
