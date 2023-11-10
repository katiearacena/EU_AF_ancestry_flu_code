#!/bin/bash

module load bedtools

#set blacklist directory
blacklist=/project2/lbarreiro/users/katie/ref/hg19_blacklist/ENCFF001TDO.bed
OUT_DIR=../H3K4me3_filtering_outputs
cat ../../H3K4me3_peaks/*.bed | sed -e 's/peak_call\///' -e 's/_H3K4me3.*//' | sort -k1,1 -k2,2n | 
bedtools merge -i - -sorted -d -174.5 -c 4,4 -o distinct,count_distinct > $OUT_DIR/cts.bed

#COUNT NUMBER OF FLU/NI
cat $OUT_DIR/cts.bed | awk -F'|' '{print gsub(/NI/,"") "\t" NR}' | awk '{print $1}' > $OUT_DIR/NI_nums.txt
cat $OUT_DIR/cts.bed | awk -F'|' '{print gsub(/Flu/,"") "\t" NR}' | awk '{print $1}' > $OUT_DIR/Flu_nums.txt
paste $OUT_DIR/cts.bed $OUT_DIR/NI_nums.txt $OUT_DIR/Flu_nums.txt > $OUT_DIR/cts_to_filter.txt

#now have a file that has number of NI and Flu pasted so can threshold on these.
#print only those lines	that are present in at least 50% of Flu	*OR* NI samples.
cat $OUT_DIR/cts_to_filter.txt | awk '{if($6 >= 14 || $7 >= 14) print $1"\t"$2"\t"$3"\t"$4"\t"$5}' |

#remove blacklisted regions.
bedtools intersect -a - -b $blacklist -v -sorted > $OUT_DIR/H3K4me3_peaks_filtered.bed

#Final check!

cat $OUT_DIR/H3K4me3_peaks_filtered.bed | awk '{if($5 <= 13) print }'
#should not print anything because of condition specific peak filtering.
