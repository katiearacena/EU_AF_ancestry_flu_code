# Epigenetic variation impacts ancestry-associated differences in the transcriptional response to influenza infection 

This repository contains the code needed to conduct main analyses for Aracena et al. 2022 ("Epigenetic variation impacts ancestry-associated differences in the transcriptional response to influenza infection")

Currently available on BioRxiv: https://www.biorxiv.org/content/10.1101/2022.05.10.491413v1

KA Aracena, YL Lin, K Luo, A Pacis, S Gona, Z Mu, V Yotova, R Sindeaux, A Pramatarova, MM Simon, X Chen, C Groza, D Lougheed, R Gregoire, D Brownlee, Y Li, X He, D Bujold, T Pastinen, G Bourque, LB Barreiro


## General dependencies
R (tested in version 3.6.3)

## Usage
1. Place the folder `Inputs`, available for download on Zenodo (XXXXX), in the desired working directory.
2. Place the folder `Main_analyses` in the same directory. 
3. Check individual dependencies at the header of each script stored in `Main_analyses` and install the necessary CRAN and Bioconductor packages.
4. Run the scripts stored in `Main_analyses` in R -- if  `Main_analyses` and `Inputs` are located in the same directory, results will populate in a folder named `Outputs`.

Note: Many scripts run R scripts for each data type and condition separately. In these cases the `.sh` script will use the `.sbatch` script to submit jobs using the `.R` script command line argument combination.

## Description of analyses 

### `1_DE_infection_modeling`
- `1_1_DE_infection`
	- Differential infection effects (Figure 1C)
- `1_2_PVE`
	- Percent variance explained by infection effects (Figure 1B)
- `1_3_GSEA`
	- Gene set enrichment analysis based on infection effects results (Figure 1D)
- `1_4_TF_activity_scores`
	- Transcription factor activity scores (Figure 1G)


### `2_ancestry_effects_modeling`
- `2_1_get_batch_age_corrected_cts`
	- Batch and age corrected counts used for ancestry and QTL mapping analyses (with the exception of methylation data)
- `2_2_PopDE`
	- Identification of features associated with genetic ancestry at both baseline and after flu infection (PopDE) (Figure 2A)
- `2_3_PopDR`
	- Identification of features in which response to flu infection is associated with genetic ancestry (PopDR)  
- `2_4_ancestry_GSEA`
	- Gene set enrichment analysis on the PopDE and PopDR results (Figure 5D)
- `2_5_ancestry_score` 
	- Ancestry scores (Figure 2C, 2D, S2C)
- `2 6_epi_priming`
	- Predicting transcriptional response using epigenetic marks at baseline (Figure 2E, S2D)

### `3_QTL_mapping`

- `3_1_SNP-QTL_mapping`
	- SNP QTL mapping (Figure 3A, S3A, S3C)
- `3_2_SNP-QTL_mapping_mash`
	- Applying mash to SNP QTL mapping results to detect condition specific effects (S3C, S3D)
- `3_3_STR-QTL_mapping`
	- STR QTL mapping (Figure 3A, S3B, S3C)
- `3_4_STR-QTL_mapping_mash`
	- Applying mash to STR QTL mapping results to detect condition specific effects (S3C, S3D)
- `3_5_SNP_STR_PVE`
	- Percent variance explained by the top SNP and STR (Figure 3B, S3D)
- `3_6_SNP-QTL_integration`
	- Detection of shared QTL (Figure 3A, 3B, 3C, 3D, S3A, S3B, S3C, S3D)
- `3_7_condition_specific_SNP-QTL_enrichment`
	- TF footprint enrichment for condition specific SNP-QTL (Figure 3E, S3E)
- `3_8_QTL_integration_enrichment`
	- TF footprint enrichment for shared QTL (Figure 4E)
- `3_9_metaQTL_analysis`
	- Gene expression and epigenetic mark levels in regions that are caQTL at baseline (Figure 4F, S4E)


### `4_cisRegression_modeling`
- ` 4_1_PopDE_cisreg`
	- PopDE effects after regressing out the top SNP and STR (used as input for `4_7_cisreg_ancestry_GSEA`)
- `4_2_PopDE_SNP-QTL`
	- Calculating predicted QTL effect sizes using top SNP
- `4_3_PopDE_STR-QTL`
	- Calculating predicted QTL effect sizes using top STR 
- `4_4_PopDE_SNP_STR-QTL`
	- Calculating predicted QTL effect sizes using both top SNP and STR (Figure 5B, S5A)
- `4_5_admixture_PVE`
	- Calculating percent variance explained by QTL on ancestry effects (Figure 5C, S5B)
- `4_6_cisreg_ancestry_score`
	- Ancestry scores after regressing out effects of top SNP and STR (Figure S5C)
- `4_7_cisreg_ancestry_GSEA`
	- Gene set enrichment analysis on the PopDE effects after regressing out the top SNP and STR (Figure 5D)
