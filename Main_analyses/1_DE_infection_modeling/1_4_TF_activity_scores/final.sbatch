#!/bin/bash

#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=54GB
#SBATCH --partition=broadwl

module load python
module load bedtools
conda activate py3.7

mkdir -p /project2/lbarreiro/users/katie/EU_AF_ancestry_flu/HINT/outputs/meta/change/change

mkdir -p $PWD/OUTPUTS

# Creates a set of footprints that overlap with the differential peaks we are subsetting on
bedtools intersect -a ALL_MERGED_FOOTPRINTS -b \
PEAKS_FILE -wa > \
$PWD/OUTPUTS/change_footprints.bed


# Does motif matching in the differential peaks
rgt-motifanalysis matching --output-location $PWD/outputs \
--organism hg19 --rand-proportion 10 --filter "species:sapiens" \
--input-files $PWD/OUTPUTS/TAG_footprints.bed

# TF activity across TFs between conditions for a single individual
rgt-hint differential --organism=hg19 --bc --nc 28 \
--mpbs-files=$PWD/OUTPUTS/TAG_footprints_mpbs.bed,$PWD/OUTPUTS/TAG_footprints_mpbs.bed \
--reads-files=NI_BAM,FLU_BAM \
--conditions=NI,Flu \
--colors=red,blue \
--output-profiles \
--output-location=$PWD/OUTPUTS/IND