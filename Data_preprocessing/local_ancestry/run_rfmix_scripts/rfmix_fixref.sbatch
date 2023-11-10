#!/bin/bash
#SBATCH --partition=lbarreiro
#SBATCH --account=pi-\lbarreiro
#SBATCH --time=6:00:00
#SBATCH --mem=48G
#SBATCH --nodes=1

cd $SLURM_SUBMIT_DIR
echo "starting job $SLURM_JOB_NAME in $SLURM_SUBMIT_DIR"

module load python
conda activate myenv

chr=`echo $SLURM_JOB_NAME | cut -d"%" -f2`
echo "$chr"

ref_dir=/project/lbarreiro/USERS/katie/ref/1000genomes_hg37_phase3/
map_dir=/project/lbarreiro/USERS/bridget/references/maps/
sample_dir=/project/lbarreiro/USERS/katie/EU_AF_ancestry_flu/local_ancestry/phasing/outputs/
out_dir=/project/lbarreiro/USERS/katie/EU_AF_ancestry_flu/local_ancestry/rfmix/outputs/

rfmix -f $sample_dir/ALL.autosomes.shapeit4.phased.1KGref.fixref.applied.vcf \
	-r $ref_dir/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
	-m $ref_dir/integrated_call_samples_v3.20130502.CEU.YRI.only.panel \
	-g $map_dir/allchrs.b37.gmap \
	-o $out_dir/ref_$chr \
	--chromosome=${chr}


## NOTE
# IF ARE CONSIDERING ALLPOPULATIONS VERSUS JUST YRI AND CEU
#-m $ref_dir/integrated_call_samples_v3.20130502.ALL.panel
