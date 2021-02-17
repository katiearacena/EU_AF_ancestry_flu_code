#!/bin/bash

############################################
#SBATCH --job-name=WGBS_MATRIXeQTL
#SBATCH --output=logs/WGBS_ME-%j.out
#SBATCH --error=logs/WGBS_ME-%j.err
#SBATCH --time=20:00:00
#SBATCH --partition=bigmem2
#SBATCH --mem=300G
#SBATCH --nodes=1
#SBATCH --mail-user=karacena@uchicago.edu
#SBATCH --mail-type=ALL
############################################

Rscript pt2.R
