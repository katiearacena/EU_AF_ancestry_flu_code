##################
# Compile HipSTR #
##################
git clone https://github.com/tfwillems/HipSTR.git
cd HipSTR
make
cd ../

#############################
# Download HipSTR reference #
#############################
wget https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/GRCh37.hipstr_reference.bed.gz

##############
# Run HipSTR #
##############
sbatch run_HipSTR.slurm

#################
# Filter output #
#################
python scripts/filter_vcf.py  --vcf                   str_calls_EU_AF.vcf.gz
                              --min-call-qual         0.9
                              --max-call-flank-indel  0.15
                              --max-call-stutter      0.15
			      --min-call-allele-bias  -2
			      --min-call-strand-bias  -2  > str_calls_EU_AF_filtered.vcf

##################################################
# Create Matrix-eQTL input from the STR vcf file #
##################################################
python create_matrixQTL_inputs_from_STRvcf.py
