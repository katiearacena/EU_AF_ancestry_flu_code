#!/bin/bash

#SBATCH --job-name=ldsc_annot_bedfiles
#SBATCH --output=ldsc_annot_bedfiles_%J.out
#SBATCH --error=ldsc_annot_bedfiles_%J.err
#SBATCH --account=pi-xinhe
#SBATCH --partition=broadwl
#SBATCH --mem=10G

## Compute binary annotations and annotation-specific LD scores using multiple bedfiles
## https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial

PREFIX_ANNOT=$1
BED_DIR=$2
PREFIX_BFILE=$3
echo "prefix of annot: ${PREFIX_ANNOT}"
echo "BED directory: ${BED_DIR}"
echo "prefix of bfile: ${PREFIX_BFILE}"

dir_LDSC=/project2/xinhe/kevinluo/ldsc/

source activate ldsc

module load R/4.2.0

mkdir -p ${dir_LDSC}/annot/ldscores/${PREFIX_ANNOT}/

cd ~/projects/flu_regulatoryQTLs/code

for chrom in {1..22}
do
  echo ${chrom}

  ## Step 1: Creating an annot file
  echo "Make ldsc-friendly annotation files for ${PREFIX_ANNOT}"
  Rscript ../code/make_ldsc_binary_annots_from_bedfiles.R \
  --beddir ${BED_DIR} \
  --bimfile ${PREFIX_BFILE}.${chrom}.bim \
  --out ${dir_LDSC}/annot/ldscores/${PREFIX_ANNOT}/${PREFIX_ANNOT}.${chrom}.annot.gz

  ## Step 2: Computing LD scores with an annot file
  echo "Computing LD scores with annot file ${PREFIX_ANNOT}.${chrom}.annot.gz"
  python $HOME/softwares/ldsc/ldsc.py \
  --l2 \
  --bfile ${PREFIX_BFILE}.${chrom} \
  --print-snps ${dir_LDSC}/LDSCORE/listHM3.txt \
  --ld-wind-cm 1 \
  --annot ${dir_LDSC}/annot/ldscores/${PREFIX_ANNOT}/${PREFIX_ANNOT}.${chrom}.annot.gz \
  --out ${dir_LDSC}/annot/ldscores/${PREFIX_ANNOT}/${PREFIX_ANNOT}.${chrom}
done

