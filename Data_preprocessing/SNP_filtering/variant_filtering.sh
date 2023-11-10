## script used to filter variants for downstream analyses

####################################
## Remove non- biallelic variants ##
####################################

module load bcftools
module load plink
module load java

DIR= <out_dir>

# Remove all sites that are not bi-allelic from input file
bcftools view -m2 -M2 -v snps <input.vcf> >> ${DIR}/bi-allelic-only.35Samples.noIPSC.hc.vqsr.vcf

#############################################
## Perform other filtering (detailed below ##
#############################################

# Using the following filters:
# 1. remove SNPs if they have a call rate of <90% across all samples (--geno 0.10)
# 2. remove SNPs if they deviate from HWE with p < 10−5 (--hwe 0.00001)
# 3. remove SNPs on with MAF (minor allele frequency) less than 5% (--maf 0.05)
# Note: because of the small sample size we do not remove any individuals (with the mind flag)
plink --vcf ${DIR}/bi-allelic-only.35Samples.noIPSC.hc.vqsr.vcf --chr 1-22 --double-id --allow-extra-chr --make-bed --out bi-allelic-only.35Samples.noIPSC.hc.vqsr

# perform filtering
plink --bfile ${DIR}/bi-allelic-only.35Samples.noIPSC.hc.vqsr --double-id --allow-extra-chr --geno 0.10 --maf 0.05 --hwe 0.00001 --make-bed --out ${DIR}/bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered

# recode into a vcf file
plink --bfile ${DIR}/bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered --recode vcf --out ${DIR}/bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered

grep -v "^#" ${DIR}/bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered.vcf | wc -l

bcftools view ${DIR}/bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered.vcf -Oz -o ${DIR}/bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered.vcf.gz
bcftools index ${DIR}/bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered.vcf.gz

## add to vcf file. ignore REF and ALT for matching
bcftools annotate -a ref/00-All.vcf.gz --collapse snps -c ID -o bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered.rs.ID.annotated.vcf.gz bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered.vcf.gz
bcftools view -Ov -o bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered.rs.ID.annotated.vcf bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered.rs.ID.annotated.vcf.gz
cut -f 1,2,3,4,5 bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered.rs.ID.annotated.vcf | sed 's/[\t]/,/g' > cols.txt

## in R
rsIDs <- read.csv("cols.txt", header=FALSE, comment.char="#")
colnames(rsIDs) <- c("chr", "pos", "rsID", "A1", "A2")
## merge chr, pos and allels into column
rsIDs$snps <- paste0(rsIDs$chr, "_", rsIDs$pos, "_", rsIDs$A1, "_", rsIDs$A2)
write.table(rsIDs, "rsIDs.for.matrixeqtl.txt")