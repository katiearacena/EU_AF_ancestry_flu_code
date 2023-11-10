#####################################
# Download codes for model training #
#####################################
git clone https://github.com/hakyimlab/PredictDBPipeline.git

##############################
# Create directory structure #
##############################
mkdir data
cd data 
mkdir output
mkdir intermediate
mkdir input
cd input
mkdir expression_phenotypes
mkdir genotypes
mkdir annotations
cd annotations
mkdir gene_annotation
mkdir snp_annotation
cd ../../../

##################################
# Parse and populate input files #
##################################
## SNP genotypes
cp ../../../Inputs/QTL_mapping/SNP_genotypes.txt data/input/genotypes/SNP_genotypes.txt

## SNP annotations
python parse_snp.py

## expression phenotypes
cp ../../../Inputs/counts_matrices/*filtered.counts.txt data/input/expression_phenotypes/

## gene(peak) annotations
cp ../../../Inputs/QTL_mapping/*positions.txt data/input/annotations/gene_annotation/
python parse_gene.py


#############################
# Modify the joblogs folder #
#############################
mv create_model.R joblogs
mv model_parameters.py joblogs 
mv job_specs.sbatch joblogs 


##############################################################
# Loop to run model training across datatypes and conditions #
##############################################################
cd joblogs
for study in RNAseq ATACseq H3K27ac H3K27me3 H3K4me1 H34Kme3
	do 
            mkdir $study
            cp * $study/
            cd $study
            for condition in NI Flu
                sed -i "s|'DATATYPE'|$study|g" model_parameters.py
                sed -i "s|'CONDITION'|$condition|g" model_parameters.py
                python preprocess.py
		do for chrom in $(seq 1 22)
			do 
				alpha=0.5
				window=1e6
				echo " *************************************************** "
				echo $study $chrom $condition
				sbatch job_specs.sbatch $study $chrom $alpha $window $condition
				echo " *************************************************** "
			done
                python post_process.py
		done
            
	done







