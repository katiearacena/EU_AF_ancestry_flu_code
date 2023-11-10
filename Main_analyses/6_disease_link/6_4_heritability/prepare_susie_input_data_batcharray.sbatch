#!/bin/bash

#SBATCH --job-name=prepare_susie_input_batches
#SBATCH --output=prepare_susie_input_batches_%A_%a.out
#SBATCH --error=prepare_susie_input_batches_%A_%a.err
#SBATCH --account=pi-xinhe
#SBATCH --partition=broadwl
#SBATCH --array=1-5
#SBATCH --mem=25G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

# --ntasks-per-node will be used in doParallel.R to specify the number
# of cores to use on the machine.

## parameters
EXP_NAME=$1

IBATCH=$(($SLURM_ARRAY_TASK_ID))

echo "Prepare SuSiE input data for: ${EXP_NAME}"

module load R/4.0.4

echo "Submit job for batch #${IBATCH}."

Rscript ~/projects/flu_regulatoryQTLs/code/prepare_susie_input_data_batches.R \
  --phenotype ${EXP_NAME} \
  --vcf /project2/xinhe/kevinluo/flu_regulatoryQTLs/genotype/202106_imputed_vcf_v2/imputed_ALL.chrs.annotated.v2.renamed.vcf.gz \
  --outdir /scratch/midway2/kaixuan/flu_regulatoryQTLs/susie/input \
  --ncores 8 \
  --nbatches 5 \
  --ibatch ${IBATCH}

echo "Done."

