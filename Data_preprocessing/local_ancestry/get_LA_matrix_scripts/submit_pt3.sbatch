#!/bin/bash
#SBATCH --partition=lbarreiro
#SBATCH --account=pi-\lbarreiro
#SBATCH --time=96:00:00
#SBATCH --mem=48G
#SBATCH --nodes=1

module load R
Rscript pt3_calc_overlap.R
