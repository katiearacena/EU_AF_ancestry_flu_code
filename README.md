# Epigenetic variation impacts individual differences in the transcriptional response to influenza infection 

This repository contains the code needed to conduct main analyses for Aracena et al., 2024.

KA Aracena, YL Lin, K Luo, A Pacis, S Gona, Z Mu, V Yotova, R Sindeaux, A Pramatarova, MM Simon, X Chen, C Groza, D Lougheed, R Gregoire, D Brownlee, C Boye, R Pique-Regi, Y Li, X He, D Bujold, T Pastinen, G Bourque, LB Barreiro Epigenetic variation impacts individual differences in the transcriptional response to influenza infection. Nat Genet 56, 408–419 (2024). https://doi.org/10.1038/s41588-024-01668-z

## General dependencies
R (tested in version 3.6.3)

## Usage
1. Place the folder `Inputs`, available for download on Zenodo (https://zenodo.org/records/10108241), in the desired working directory.
2. Place the folders `Data_Preprocessing` and `Main_analyses` in the same directory. 
3. Check individual dependencies at the header of each script and install the necessary CRAN and Bioconductor packages.
4. Run the scripts stored in `Main_analyses` in R -- if  `Main_analyses` and `Inputs` are located in the same directory, results will populate in a folder named `Outputs`.

Note: Many scripts run R scripts for each data type and condition separately. In these cases the `.sh` script will use the `.sbatch` script to submit jobs using the `.R` script command line argument combination.

## Description of `Main_analyses` directories 

### `1_DE_infection_modeling`
- `1_1_DE_infection`
	- Differential infection effects (Figure 1C)
- `1_2_PVE`
	- Percent variance explained by infection effects (Figure 1B)
- `1_3_GSEA`
	- Gene set enrichment analysis based on infection effects results (Figure 1D)
- `1_4_TF_activity_scores`
	- Transcription factor activity scores (Figure 1G)
- `1_5_PVE_mock_comparison`
	-Percent variance explained by technical variation (mock v NI samples; Extended Data Fig. 1A)
- `1_6_mock_correlation`
	- Correlation of read counts for mock v. NI samples (Supplementary Figure 2B)



### `2_ancestry_effects_modeling`
- `2_1_get_batch_age_corrected_cts`
	- Batch and age corrected counts used for ancestry and QTL mapping analyses (with the exception of methylation data)
- `2_2_PopDE`
	- Identification of features associated with genetic ancestry at both baseline and after flu infection (PopDE) (Figure 2A)
- `2_3_PopDR`
	- Identification of features in which response to flu infection is associated with genetic ancestry (PopDR)  
- `2_4_ancestry_GSEA`
	- Gene set enrichment analysis on the PopDE and PopDR results (Figure 6D)
- `2_5_ancestry_score` 
	- Ancestry scores (Figure 2C, 2D)
- `2 6_epi_priming`
	- Predicting transcriptional response using epigenetic marks at baseline (Figure 2E, Extended Data Fig. 2F)
- `2_7_PopDE_local_ancestry`
	-PopDE effects using local ancestry estimates (Extended Data Fig. 2A)
- `2_8_ancestry_scores_FDR_changed`
	-Show ancestry scores are robust to various FDR thresholds (Extended Data Fig. 2C-E)


### `3_QTL_mapping`

- `3_1_SNP-QTL_mapping`
	- SNP QTL mapping (Figure 3A, Extended Data Fig. 4A)
- `3_2_SNP-QTL_mapping_mash`
	- Applying mash to SNP QTL mapping results to detect condition specific effects (Extended Data Fig. 4C-D)
- `3_3_STR-QTL_mapping`
	- STR QTL mapping (Figure 3A, Extended Data Fig. 4B)
- `3_4_STR-QTL_mapping_mash`
	- Applying mash to STR QTL mapping results to detect condition specific effects (Extended Data Fig. 4C-D)
- `3_5_SNP_STR_PVE`
	- Percent variance explained by the top SNP and STR (Figure 3B, Extended Data Fig. 4D)
- `3_6_SNP-QTL_integration`
	- Detection of shared QTL (Figure 4A-D, Extended Data Fig. 5A-D)
- `3_7_condition_specific_SNP-QTL_enrichment`
	- TF footprint enrichment for condition specific SNP-QTL (Figure 3E, Extended Data Fig. 4E)
- `3_8_QTL_integration_enrichment`
	- TF footprint enrichment for shared QTL (Figure 4E)
- `3_9_metaQTL_analysis`
	- Gene expression and epigenetic mark levels in regions that are caQTL at baseline (Figure 5A-C, Extended Data Fig. 6A-B)
- `3_9_REMAPPING_QTL_ANALYSIS`
	-QTL mapping using re-mapped read counts (Supplementary Figure 3D-E)


### `4_cisRegression_modeling`
- ` 4_1_PopDE_cisreg`
	- PopDE effects after regressing out the top SNP and STR (used as input for `4_7_cisreg_ancestry_GSEA`)
- `4_2_PopDE_SNP-QTL`
	- Calculating predicted QTL effect sizes using top SNP
- `4_3_PopDE_STR-QTL`
	- Calculating predicted QTL effect sizes using top STR 
- `4_4_PopDE_SNP_STR-QTL`
	- Calculating predicted QTL effect sizes using both top SNP and STR (Figure 6B, Extended Data Fig. 7A)
- `4_5_admixture_PVE`
	- Calculating percent variance explained by QTL on ancestry effects
- `4_6_cisreg_ancestry_score`
	- Ancestry scores after regressing out effects of top SNP and STR (Extended Data Fig. 7B)
- `4_7_cisreg_ancestry_GSEA`
	- Gene set enrichment analysis on the PopDE effects after regressing out the top SNP and STR (Figure 6C)


### `5_enhancerRNAs`
- ` 5_1_DE_infection_modeling`
	- Infection effects for identified enhancerRNAs
- ` 5_2_popDE_modeling`
	- PopDE modeling for enhancerRNAs
- ` 5_3_eRNA_QTL_mapping`
	-SNP QTL mapping with enhancerRNAs
- ` 5_4_eRNA_QTL_integration`
	-Identifying enhancerRNA QTL and sharing with other data types 
- ` 5_5_eRNA_mash`
	-Applying mash to QTL effects


### `6_disease_link`
- ` 6_1_Predixcan_model_training`
	- Training models for Predixcan analysis 
- ` 6_2_Summary-Predixcan_of_GWAS`
	Running Predixcan for GWAS (Figure 7E, Extended Data Fig. 8B)
- ` 6_3_colocalization`
	Scripts to run colocalization analysis (Figure 7A-B, Extended Data Fig. 8A)
- ` 6_4_heritability`
	Scripts to run heritability analyses (Figure 7C-D, Extended Data Fig. 9A-D, Supplementary Figure 4B)

