## Load libraries
library(plyr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(reshape2)
library(tidyr)
library(limma)
library(DSS)

#################################
##### LOAD COMMAND LINE ARGS ####
#################################

args<-commandArgs(TRUE)
#datatype (for directory structure)
data <- args[1]

if(length(args)==0)
{
  print("WARNING: No arguments supplied.")
}
print(args)

##############################################
## Create directory structure and load data ##
##############################################

folder = "4_cisregression_modeling"
## Set directory 3 steps above script.
setwd('../../../')

## Create directory structure to save outputs.
system(paste0("mkdir -p Outputs/",folder,"/predicted_pop_diffs_STRs/", data,"/"))
out_dir <- paste0("Outputs/",folder,"/predicted_pop_diffs_STRs/", data,"/")

## read in meta data
meta_data <- read.table(paste0("Inputs/metadata/", data, "_metadata.txt"))
meta_data$Genotyping_ID <- as.factor(meta_data$Genotyping_ID)
meta_data$Condition = factor(meta_data$Condition, levels=c("NI","Flu"))

## Read in popDE results 
popDE_NI_genes <- read.table(paste0("Outputs/2_ancestry_effects_modeling/PopDE_analysis/", data, "/results/popDE_NI.txt"))
popDE_FLU_genes <- read.table(paste0("Outputs/2_ancestry_effects_modeling/PopDE_analysis/", data, "/results/popDE_FLU.txt"))

## Subset on significant ones (10% qval)
popDE_NI_genes_sig <- popDE_NI_genes[popDE_NI_genes$qvals_NI < .10, ]
popDE_FLU_genes_sig <- popDE_FLU_genes[popDE_FLU_genes$qvals_Flu < .10, ]


## read in STR genotypes
STRgenotypes = read.table(paste0("Inputs/QTL_mapping/STR_genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)
STRgenotypes$str <- rownames(STRgenotypes)

## load matrixeqtl best strs files for each condition
eSTR_NI_results <- read.table(paste0("Outputs/3_QTL_mapping/STR-QTL_mapping/",data,"/NI/", data,"_NI/results_best_STRs_with_qval.txt"))
eSTR_FLU_results <- read.table(paste0("Outputs/3_QTL_mapping/STR-QTL_mapping/",data,"/Flu/", data,"_Flu/results_best_STRs_with_qval.txt"))

STR_FLU_beta_SE_matrix_ALL <- arrange(eSTR_FLU_results, gene)
STR_NI_beta_SE_matrix_ALL <- arrange(eSTR_NI_results, gene)

## Will be true if both dfs are in the same order.
all(STR_FLU_beta_SE_matrix_ALL$gene == STR_NI_beta_SE_matrix_ALL$gene)

## Bind and rename column ids
STR_all_conditions <- cbind(STR_FLU_beta_SE_matrix_ALL, STR_NI_beta_SE_matrix_ALL)
colnames(STR_all_conditions) <- c("strs_flu", "gene_flu", "statistic_flu", "pvalue_flu", "FDR_flu", "flu_beta", "flu_qval", 
                                  "strs_NI", "gene_NI", "statistic_NI", "pvalue_NI", "FDR_NI", "NI_beta", "NI_qval")

eSTR_Results <- select(STR_all_conditions, "gene_flu", "strs_flu", "strs_NI", "flu_beta", "NI_beta")
colnames(eSTR_Results) <- c("gene", "Flu_strs", "NI_strs", "Flu_beta", "NI_beta")



###############################################################
## Loop to calculate predicted using QTL effect and genotype ##
###############################################################

infection <- c("NI", "Flu")

for (i in 1:length(infection)){
  
  infection_i <- infection[i]
  

  ## subset matrix eSTR results
  eSTR_Results_subset <- select(eSTR_Results, paste0(infection_i, "_beta"), paste0("gene"), paste0(infection_i, "_strs"))
  colnames(eSTR_Results_subset)<- c("beta", "gene", "strs")
  
  ## subset eQTL results on top STRs and popDE genes
  ## this will exclude any sig popDE features that do not have a QTL and those that were not included in QTL testing
  ## (those that are on X, Y, MT contigs for all datatypes and cpg sites that or have low coverage and 0 variation for WGBS)
  if (infection_i=="NI"){
    eSTR_Results_subset <- eSTR_Results_subset[eSTR_Results_subset$gene %in% rownames(popDE_NI_genes_sig), ]
  } else {
    eSTR_Results_subset <- eSTR_Results_subset[eSTR_Results_subset$gene %in% rownames(popDE_FLU_genes_sig), ]
  }
  
  ## pull the genotypes for these SNPgenes across individuals
  STRgenotypes_subset <- STRgenotypes[rownames(STRgenotypes) %in% eSTR_Results_subset$strs, ]
  eSTR_Results_subset <- eSTR_Results_subset[!duplicated(eSTR_Results_subset$strs),]
  dim(STRgenotypes_subset)[1]==dim(eSTR_Results_subset)[1]
  
  ## change -9 to NA
  STRgenotypes_subset[STRgenotypes_subset == -9] <- NA
  
  ## put genotypes and eQTL results in order
  STRorder_vec <- as.character(eSTR_Results_subset$strs)
  STRgenotypes_subset <- STRgenotypes_subset[match(STRorder_vec, rownames(STRgenotypes_subset)),]
  STRgenotypes_subset$str <- NULL
  ## are in same order if = 0.
  length(which(eSTR_Results_subset$snps!=rownames(STRgenotypes_subset)))
  
  ## multiply genotypes by effect sizes
  ## for SNPs: genotype (0,1, or 2) * beta
  STRpopValues <- sweep(STRgenotypes_subset, MARGIN = 1, eSTR_Results_subset$beta, `*`)
  
  
  
  
  
  ## create predicted exp matrix
  predicted_expression_matrix <- STRpopValues
  predicted_expression_matrix$strs <- rownames(predicted_expression_matrix)
  predicted_expression_matrix <- join(predicted_expression_matrix, eSTR_Results_subset, by = "strs")
  rownames(predicted_expression_matrix) <- predicted_expression_matrix$gene
  predicted_expression_matrix$gene <- NULL
  predicted_expression_matrix$strs <- NULL
  predicted_expression_matrix$beta <- NULL
  colnames(predicted_expression_matrix) <- paste0(colnames(predicted_expression_matrix),"_",infection_i)
  
  ## write these as the predicted expression matrices
  write.table(predicted_expression_matrix, paste0(out_dir,infection_i,"_predicted_expression.txt"), quote = FALSE)
  
  ###########################
  ##### RUN POPDE MODEL #####
  ###########################
  
  #subset to only include those samples for which there is a 1 in the "PopDEset" column
  meta_data=meta_data[which(meta_data$PopDE_set==1),]
  #ensure that the Sample_ID are the rownames of the cols dataframe
  rownames(meta_data)=meta_data$Sample
  #ensure that the cols dataframe is ordered alphabetically (based on sample_id)
  meta_data=meta_data[order(rownames(meta_data)),]
  
  rownames(meta_data)<- sub("_H3K27ac", "", rownames(meta_data))
  rownames(meta_data)<- sub("_H3K4me1", "", rownames(meta_data))
  rownames(meta_data)<- sub("_H3K27me3", "", rownames(meta_data))
  rownames(meta_data)<- sub("_H3K4me3", "", rownames(meta_data))
  
  cols <- meta_data[rownames(meta_data) %in% colnames(predicted_expression_matrix),]
  predicted_expression_matrix <- predicted_expression_matrix[, colnames(predicted_expression_matrix) %in% rownames(cols)]
  
  length(which(colnames(predicted_expression_matrix)!=rownames(cols)))
  
  ## run model for other datatypes
  design = model.matrix(~Admixture, data = cols)
  fit <-lmFit(predicted_expression_matrix,design)
  fit <- eBayes(fit)
  
  betas = as.data.frame(fit$coefficients[, which(colnames(fit$coefficients) %in% c("Admixture"))]); colnames(betas)[1] <- "betas"
  
  ## write outs
  write.table(betas, paste0(out_dir,infection_i,"_predicted_PopDE_betas.txt"), quote = FALSE)
}

