OUT_DIR=../../H3K4me1_filtering_outputs
ALIGNMENTS=/project2/lbarreiro/users/katie/EU_AF_ancestry_flu/data/H3K4me1_alignment/not_used_inChromHMM/not_used_in_any_analysis

#submit jobs
suffix='sorted.dup.bam'
for f in `cat $OUT_DIR/mock_samples.txt`; do f=${f/\.$suffix/};
echo '#!/bin/sh'> $f.sh;
echo "module load gcc/7.2.0
/project2/lbarreiro/programs/subread-1.6.4-Linux-x86_64/bin/featureCounts -T 12 -p -P -d 30 -D 1000 --ignoreDup -a $OUT_DIR/H3K4me1_peaks_filtered.mod.gtf $ALIGNMENTS/$f.$suffix -o $OUT_DIR/$f.counts.paired.noDup.txt" >> $f.sh ; sbatch -p broadwl -D $PWD -o $f.out -J $f --time=10:00 --mem=48G -N 1 -n 12 $f.sh; done


#combine all files into 1 count file.

cd $OUT_DIR

mv *.noDup.txt *.txt.summary raw_counts
mv ../scripts/GET_MOCK/*.out out_logs


cd raw_counts
paste -d "\t" $(ls -1v *.counts.paired.noDup.txt) | sed '/^#/d' | cut -f 1,`seq --separator="," 7 7 462` | sed -e 's/.sorted.dup.bam//g' -e 's/Geneid/PeakID/g' -e 's/alignment\///g' > H3K4me1_peaks_with_MOCK.filtered.paired.noDup.counts


#to remove some of the extra id label information in the counts file.
sed -i "s|/project2/lbarreiro/users/katie/EU_AF_ancestry_flu/data/||g" H3K4me1_peaks_with_MOCK.filtered.paired.noDup.counts
sed -i "s|H3K4me1_||g" H3K4me1_peaks_with_MOCK.filtered.paired.noDup.counts

mv H3K4me1_peaks_with_MOCK.filtered.paired.noDup.counts /project2/lbarreiro/users/katie/EU_AF_ancestry_flu/mock_counts


