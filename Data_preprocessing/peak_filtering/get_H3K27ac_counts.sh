#!/bin/bash

## Example code to get samples text file
##ls *.bed | sed -e 's/_peaks.narrowPeak.bed//' > ../FILTER_PEAKS/H3K4me3_filtering_outputs/samples.txt


##run in scripts directory##

OUT_DIR=../H3K27ac_filtering_outputs
ALIGNMENTS=/project2/lbarreiro/users/katie/EU_AF_ancestry_flu/data/H3K27ac_alignment

python bed2gtf.py -i $OUT_DIR/H3K27ac_peaks_filtered.bed -o $OUT_DIR/H3K27ac_peaks_filtered.gtf 
sed -i -e 's/peak_id/gene_id/' -e 's/peak/exon/' $OUT_DIR/H3K27ac_peaks_filtered.gtf

#remove any rows that have a start position=0 (if you do not this will throw an error for feature counts).
awk '($4>0)' $OUT_DIR/H3K27ac_peaks_filtered.gtf > $OUT_DIR/H3K27ac_peaks_filtered.mod.gtf

#submit jobs
suffix='sorted.dup.bam'
for f in `cat $OUT_DIR/samples.txt`; do f=${f/\.$suffix/}; 
echo '#!/bin/sh'> $f.sh;
echo "module load gcc/7.2.0
/project2/lbarreiro/programs/subread-1.6.4-Linux-x86_64/bin/featureCounts -T 12 -p -P -d 30 -D 1000 --ignoreDup -a $OUT_DIR/H3K27ac_peaks_filtered.mod.gtf $ALIGNMENTS/$f.$suffix -o $OUT_DIR/$f.counts.paired.noDup.txt" >> $f.sh ; sbatch -p broadwl -D $PWD -o $f.out -J $f --time=10:00 --mem=48G -N 1 -n 12 $f.sh; done


#combine all files into 1 count file.

cd $OUT_DIR
paste -d "\t" $(ls -1v *.counts.paired.noDup.txt) | sed '/^#/d' | cut -f 1,`seq --separator="," 7 7 462` | sed -e 's/.sorted.dup.bam//g' -e 's/Geneid/PeakID/g' -e 's/alignment\///g' > H3K27ac_peaks.filtered.paired.noDup.counts
#store raw counts
mkdir raw_counts
mkdir out_logs
mv *.noDup.txt *.txt.summary raw_counts
mv ../scripts/*.out out_logs

#to remove some of the extra id label information in the counts file.
sed -i "s|/project2/lbarreiro/users/katie/EU_AF_ancestry_flu/data/||g" H3K27ac_peaks.filtered.paired.noDup.counts
sed -i "s|_H3K27ac||g" H3K27ac_peaks.filtered.paired.noDup.counts
