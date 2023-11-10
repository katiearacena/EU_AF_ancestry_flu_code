## command to submit script
for c in {1..2}; do sbatch --job-name rmfix%$c rfmix_fixref.sbatch; done



## TO SUBSET BASED ON POPULATIONS. CREATE REF VCF WITH APPROPRIATE SUBPOPS INCLUDED AND OTHERS EXCLUDED
## PERFORMED IN REF_DIR
## add # to all rows/lines
sed 's/^/#/' integrated_call_samples_v3.20130502.ALL.panel > integrated_call_samples_v3.20130502.subset.panel
## remove # from CEU or YRI lines only 
awk '{if ($2=="YRI" || $2=="CEU") { gsub("#","")} else {}}1' integrated_call_samples_v3.20130502.subset.panel > integrated_call_samples_v3.20130502.CEU.YRI.only.panel
## make sure only tabs, no spaces
#awk -v OFS="\t" '$1=$1' integrated_call_samples_v3.20130502.CEU.YRI.only.panel

wc integrated_call_samples_v3.20130502.ALL.panel
wc integrated_call_samples_v3.20130502.CEU.YRI.only.panel


