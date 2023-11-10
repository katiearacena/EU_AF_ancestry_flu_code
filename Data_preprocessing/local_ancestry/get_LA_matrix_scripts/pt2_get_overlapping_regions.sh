## Written by KA on 1/12/23
## script is part 2 to get RFmix LA results into the correct format for popDE LA analysis
## the inputs are as follows:

## 1. LA_matrix_for_popDE.txt which has collpased LA estimates for each individual (cols) for each LA region (rows). 
## this matrix has removed haplotype info since it's no longer relevant so only options are 0,1 and 2 for each individual
## 0 - 2 CEU copies, 2- 2 YRI copies 

## 2. positions file (${DAT}_positions.txt). I used one I create in QTL mapping file. In this file col1=region ID, col2=chromosome, col3=startpos, col4=endpos
## (make sure "chr" or no chr is consistent between files)

## enter working directory
cd /project/lbarreiro/USERS/katie/EU_AF_ancestry_flu/local_ancestry/popDE_inputs

# not needed on midway3 
# module load bedtools

DATA_TYPES=(RNAseq ATACseq H3K27ac H3K27me3 H3K4me1 H3K4me3 WGBS)

## extract appropriate columns for LA regions
awk -v OFS='\t' '{print $2,$3,$4}' LA_matrix_for_popDE.txt | sed '1d' > LA_regions.bed
## sort bed file
sort -k1,1 -k2,2n LA_regions.bed > LA_regions.sorted.bed

##### PREPARE FILE FOR EACH PHENOTYPE

LEN=${#DATA_TYPES[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       DAT=${DATA_TYPES[$NUM]}
       echo "running for "${DAT}

       ## first convert positions file into initial bed
       awk '{print $3,$4,$5}' ${DAT}_positions.txt | sed '1d' > ${DAT}_pos.bed
       ## format bed file correctly
       ## remove "chr", quotes and add tab delimination
       sed -i 's/chr//g' ${DAT}_pos.bed
       sed -i 's/"//g' ${DAT}_pos.bed
       sed -i 's/ /\t/g' ${DAT}_pos.bed

       ## sort bed file
       sort -k1,1 -k2,2n ${DAT}_pos.bed > ${DAT}_positions.sorted.bed
       ##delete unsorted
       rm ${DAT}_pos.bed

       ## add +/100 KB to phenotype positions file.
       awk -v s=100000 '{print $1"\t"$2-s"\t"$3+s"\t"$1"\t"$2"\t"$3}' ${DAT}_positions.sorted.bed > ${DAT}_positions.sorted.100kb_window.bed
       ## force all negative numbers to 0
       awk -v OFS='\t' '$2<0 {$2=0} 1' ${DAT}_positions.sorted.100kb_window.bed > ${DAT}_positions.sorted.100kb_window.nonegs.bed

       bedtools intersect -a ${DAT}_positions.sorted.100kb_window.nonegs.bed -b LA_regions.sorted.bed -wo > ${DAT}_intersected_LA_coordinates.bed

done

## now have file that has overlap for all features with LA.
## cols 1-3 are coordinate of 100+/ kb region, cols 4-6 are original feature region, cols 7-9 are overlapping LA regions

