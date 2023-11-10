# command to submit phasing script
for c in {1..2}; do sbatch --job-name bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered.rs.ID.annotated.v4.vcf.gz%chr$c shapeit4.sbatch; done

## to combine all chrs
bcftools concat bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered.rs.ID.annotated.v4.chr{1..22}.shapeit4.phased.1KGref.vcf.gz -Oz -o  ALL.autosomes.shapeit4.phased.1KGref.vcf.gz
bcftools concat bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered.rs.ID.annotated.v4.chr{1..22}.shapeit4.phased.NOref.vcf.gz -Oz -o  ALL.autosomes.shapeit4.phased.NOref.vcf.gz

