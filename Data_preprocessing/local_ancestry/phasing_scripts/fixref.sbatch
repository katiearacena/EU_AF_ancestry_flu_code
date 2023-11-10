#!/bin/bash
#SBATCH --partition=lbarreiro
#SBATCH --account=pi-\lbarreiro
#SBATCH --time=96:00:00
#SBATCH --mem=48G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20

### on midway3
module load htslib
module load java
module load plink

DIR=/project/lbarreiro/USERS/katie/EU_AF_ancestry_flu/local_ancestry/phasing/scripts
FILE=bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered.rs.ID.annotated.v4
REF=/project/lbarreiro/USERS/katie/ref/hg37_genpipes/Homo_sapiens.GRCh37.fa

ONEG=/project/lbarreiro/USERS/katie/ref/1000genomes_hg37_phase3/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes

# check if there are mismatched samples
bcftools +fixref ${DIR}/${FILE}.vcf -- -f ${REF}
# 20.7% mismatched
# this is due to plink's filtering which automatically marks the minor allele in the vcf as Alt allele

## now perform actual flipping
bcftools +fixref ${DIR}/${FILE}.vcf -Ov -o ${FILE}-post_fixref.vcf -- -d -f ${REF}

## check again to be sure they match
bcftools +fixref ${DIR}/${FILE}-post_fixref.vcf -- -f ${REF}
## now 100% match and 0% mismatch
 
## now for 1000gs samples
bcftools +fixref ${ONEG}.vcf.gz -Ov -o ${ONEG}-post_fixref.vcf.gz -- -d -f ${REF}
## check again to be sure they match
bcftools +fixref ${ONEG}-post_fixref.vcf.gz -- -f ${REF}
## 100% matched so there is no post-fix ref file created. 

