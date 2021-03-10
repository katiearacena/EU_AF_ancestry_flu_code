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

## Set functions
# This function mean centers the data.
mean_center_clean=function(cols_local,column){
  index=which(colnames(cols_local)==column)
  cols_local[,index]=cols_local[,index]-mean(cols_local[,index])
  return(cols_local)
}
# This function applies mean centers function to desired metadata column. Here we apply it to age.
pretty_up_cols=function(cols){

  cols=refactorize(cols)
  cols=mean_center_clean(cols,"Age")
  return(cols)
}
# This function makes sure that desired metadata columns are factors.
refactorize=function(cols_local){

  cols_local$Genotyping_ID=factor(cols_local$Genotyping_ID,levels=unique(as.character(cols_local$Genotyping_ID)))
  cols_local$Sample=factor(cols_local$Sample)
  cols_local$Condition=factor(cols_local$Condition)
  cols_local$Batch=factor(cols_local$Batch)
  return(cols_local)
}

##############################################
## Create directory structure and load data ##
##############################################

folder = "4_cisregression_modeling"
## Set directory 3 steps above script.
setwd('../../../')

## Create directory structure to save outputs.
system(paste0("mkdir -p Outputs/",folder,"/predicted_pop_diffs_SNPs", data,"/"))
out_dir <- paste0("Outputs/",folder,"/predicted_pop_diffs_SNPs/", data,"/")

## read in meta data
meta_data <- read.table(paste0("Inputs/metadata/", data, "_metadata.txt"))
meta_data$Genotyping_ID <- as.factor(meta_data$Genotyping_ID)
meta_data$Condition = factor(meta_data$Condition, levels=c("NI","flu"))

## Read in popDE results 
popDE_NI_genes <- read.table(paste0("Outputs/2_ancestry_effects_modeling/PopDE_analysis/", data, "/results/popDE_NI.txt"))
popDE_FLU_genes <- read.table(paste0("Outputs/2_ancestry_effects_modeling/PopDE_analysis/", data, "/results/popDE_FLU.txt"))

## Subset on significant ones (10% qval)
popDE_NI_genes_sig <- popDE_NI_genes[popDE_NI_genes$qvals_NI < .10, ]
popDE_FLU_genes_sig <- popDE_FLU_genes[popDE_FLU_genes$qvals_Flu < .10, ]

## read in genotypes
genotypes = read.table(paste0("Inputs/QTL_mapping/SNP_genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)
genotypes$snp <- rownames(genotypes)

## read in mash QTL results
eQTL_posterior_means <- read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_mash/", data, "/posteriorMeans.txt"), header = TRUE)
## only need results from one condition to get gene and snp info
results_both <- read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_mash/", data, "/results_BOTH.mash.txt"), header = TRUE)
colnames(eQTL_posterior_means) <- c("pm_Flu", "pm_NI")

## get column that is shared in both
eQTL_posterior_means$feature_id <- rownames(eQTL_posterior_means) 
results_both$feature_id <- rownames(results_both) 
## merge
eQTL_Results <- full_join(results_both, eQTL_posterior_means, by="feature_id")

## keep only snps, gene and pms
eQTL_Results <- select(eQTL_Results, "gene_flu", "snps_flu", "pm_Flu", "pm_NI")
colnames(eQTL_Results) <- c("gene", "snps", "pm_Flu", "pm_NI")

# ## load matrixeqtl betas
# eQTL_NI_results <- read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/NI/", data, "_NI/raw_results/result_original.txt"),header=TRUE)
# eQTL_FLU_results <- read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/Flu/", data, "_Flu/raw_results/result_original.txt"),header=TRUE)
# results_both <- full_join(eQTL_NI_results, eQTL_FLU_results, by="gene")
# results_both <- select(results_both, gene, qvalues.x, qvalues.y, snps.x, snps.y)
# colnames(results_both) <- c("gene", "qval_NI", "qval_Flu", "snps_NI", "snps_Flu")
# rownames(results_both) <- results_both$gene


# colnames(eQTL_posterior_means) <- c("pm_Flu", "pm_NI")
# eQTL_Results <- merge(results_both, eQTL_posterior_means, by=0)
# ## keep only snps, gene and pms
# eQTL_Results <- select(eQTL_Results, "gene", "snps_NI", "snps_Flu", "pm_Flu", "pm_NI", "qval_NI", "qval_Flu")

###############################################################
## Loop to calculate predicted using QTL effect and genotype ##
###############################################################

infection <- c("NI", "Flu")

for (i in 1:length(infection)){

  infection_i <- infection[i]

  ## subset matrix eQTL results
  eQTL_Results_subset <- select(eQTL_Results, paste0("pm_", infection_i), paste0("gene"), paste0("snps"))
  colnames(eQTL_Results_subset)<- c("beta", "gene", "snps")

  ## subset eQTL results on top SNPs and popDE genes
  ## this will exclude any sig popDE features that do not have a QTL and those that were not included in QTL testing
  ## (those that are on X, Y, MT contigs for all datatypes and cpg sites that or have low coverage and 0 variation for WGBS)
  if (infection_i=="NI"){
    eQTL_Results_subset <- eQTL_Results_subset[eQTL_Results_subset$gene %in% rownames(popDE_NI_genes_sig), ]
  } else {
    eQTL_Results_subset <- eQTL_Results_subset[eQTL_Results_subset$gene %in% rownames(popDE_FLU_genes_sig), ]
  }

  ## pull the genotypes for these SNPgenes across individuals
  genotypes_subset <- genotypes[rownames(genotypes) %in% eQTL_Results_subset$snps, ]
  eQTL_Results_subset <- eQTL_Results_subset[!duplicated(eQTL_Results_subset$snps),]
  dim(genotypes_subset)[1]==dim(eQTL_Results_subset)[1]

  ## change -9 to NA
  genotypes_subset[genotypes_subset == -9] <- NA

  ## put genotypes and eQTL results in order
  order_vec <- as.character(eQTL_Results_subset$snps)
  genotypes_subset <- genotypes_subset[match(order_vec, rownames(genotypes_subset)),]
  genotypes_subset$snp <- NULL
  ## are in same order if = 0.
  length(which(eQTL_Results_subset$snps!=rownames(genotypes_subset)))

  ## multiply genotypes by effect sizes
  ## for SNPs: genotype (0,1, or 2) * beta
  popValues <- sweep(genotypes_subset, MARGIN = 1, eQTL_Results_subset$beta, `*`)

  ## create predicted exp matrix
  predicted_expression_matrix <- popValues
  predicted_expression_matrix$snps <- rownames(predicted_expression_matrix)
  predicted_expression_matrix <- join(predicted_expression_matrix, eQTL_Results_subset, by = "snps")
  rownames(predicted_expression_matrix) <- predicted_expression_matrix$gene
  predicted_expression_matrix$gene <- NULL
  predicted_expression_matrix$snps <- NULL
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

if(data == "WGBS"){

    ## Mean center age & refactorize 
    cols=pretty_up_cols(cols)

    ## run the linear model 
    design = model.matrix(~Admixture + Age + Batch, data = cols)
    fit <-lmFit(predicted_expression_matrix,design)
    fit <- eBayes(fit)

    betas = as.data.frame(fit$coefficients[, which(colnames(fit$coefficients) %in% c("Admixture"))]); colnames(betas)[1] <- "betas"

} else {
    
    ## run model for other datatypes
    design = model.matrix(~Admixture, data = cols)
    fit <-lmFit(predicted_expression_matrix,design)
    fit <- eBayes(fit)

    betas = as.data.frame(fit$coefficients[, which(colnames(fit$coefficients) %in% c("Admixture"))]); colnames(betas)[1] <- "betas"
    } 

  ## write outs
  write.table(betas, paste0(out_dir,infection_i,"_predicted_PopDE_betas.txt"), quote = FALSE)
}

