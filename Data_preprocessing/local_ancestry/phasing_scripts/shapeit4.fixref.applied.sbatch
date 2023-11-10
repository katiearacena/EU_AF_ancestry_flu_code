#!/bin/bash
#SBATCH --partition=lbarreiro
#SBATCH --account=pi-\lbarreiro
#SBATCH --time=2:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20

cd $SLURM_SUBMIT_DIR

module load htslib

echo "starting job $SLURM_JOB_NAME in $SLURM_SUBMIT_DIR"

vcf=`echo $SLURM_JOB_NAME | cut -d"%" -f1`
chr=`echo $SLURM_JOB_NAME | cut -d"%" -f2`

reg=${chr/chr/}
mapdir=/project/lbarreiro/USERS/bridget/references/maps/impute2-1000GP_Phase3 
refdir=/project/lbarreiro/USERS/katie/ref/1000genomes_hg37_phase3  
pre=${vcf/\.vcf\.gz/}
shapeit_dir=/project/lbarreiro/SHARED/PROGRAMS/shapeit4-4.2.2/bin/


${shapeit_dir}\shapeit4.2 \
--input $vcf \
--region $reg \
-M $mapdir/genetic_map_${chr}_combined_b37.txt \
--reference ${refdir}/ALL.${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
--thread 20 \
-O ${pre}.${chr}.shapeit4.phased.1KGref.fixref.applied.vcf.gz

