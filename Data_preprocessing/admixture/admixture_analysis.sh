#!/bin/bash

## Script used to perform admixture

module load bcftools
module load plink
module load java

######################################
## Use harmonized vcf files as input #
#####################################

DIR= <directory>
FILE= <harmonized_vcf_file>

#Merge samples file with reference panel of choice (YRI and CEU).
bcftools merge  ${DIR}/${FILE}.vcf.bgz ${DIR}/CEU.vcf.bgz ${DIR}/YRI.vcf.bgz  -Oz -o ${DIR}/samples_with_1000genomes_post_harmonization.vcf.bgz

####################################
########### LD pruning #############
####################################

plink --vcf ${DIR}/samples_with_1000genomes_post_harmonization.vcf.bgz --double-id  --allow-extra-chr --make-bed --out ${DIR}/samples_with_1000genomes_post_harmonization
plink --bfile ${DIR}/samples_with_1000genomes_post_harmonization --allow-extra-chr --indep-pairwise 50 10 0.1
plink --bfile ${DIR}/samples_with_1000genomes_post_harmonization --allow-extra-chr --extract plink.prune.in --make-bed --out ${DIR}/samples_with_1000genomes_post_harmonization-pruned

#remove all non-integers including "chr" and the scaffold ids using the following
sed 's/^chrM\s/25\t/g; s/^chrX\s/23\t/g; s/^chrY\s/24\t/g; s/^chr//g' ${DIR}/samples_with_1000genomes_post_harmonization-pruned.bim > ${DIR}/fixed-samples_with_1000genomes_post_harmonization-pruned.bim
sed -e "s/^GL//" ${DIR}/fixed-samples_with_1000genomes_post_harmonization-pruned.bim > ${DIR}/fixed1-samples_with_1000genomes_post_harmonization-pruned.bim
sed 's/\.1//' ${DIR}/fixed1-samples_with_1000genomes_post_harmonization-pruned.bim > ${DIR}/fixed2-samples_with_1000genomes_post_harmonization-pruned.bim
cp ${DIR}/fixed2-samples_with_1000genomes_post_harmonization-pruned.bim  ${DIR}/samples_with_1000genomes_post_harmonization-pruned.bim

####################################
########## Run ADMIXTURE ###########
####################################

#Now it is time to run supervised admixture
#You must create a .pop file to run admixture which must match the other file names.
#Create samples_with_1000genomes_post_harmonization-pruned.pop

#RUN ADMIXTURE
/project2/lbarreiro/programs/admixture_linux-1.3.0/admixture --supervised --cv -j8 samples_with_1000genomes_post_harmonization-pruned.bed 2

#RUN PCA
plink --bfile samples_with_1000genomes_post_harmonization-pruned --allow-extra-chr --pca --out samples_with_1000genomes_post_harmonization-pruned-PCA

#Run HW
plink --bfile samples_with_1000genomes_post_harmonization-pruned  --hardy --allow-extra-chr --out samples_with_1000genomes_post_harmonization-pruned-HW
grep ALL samples_with_1000genomes_post_harmonization-pruned-HW.hwe > samples_with_1000genomes_post_harmonization-pruned-HW.hwe_reduced