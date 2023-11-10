# Load libraries
library(limma)
library(edgeR)
library(ggplot2)
library(cowplot)
library(reticulate)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(gridExtra)
library(ggfortify)
library(qvalue)

#################################
##### LOAD COMMAND LINE ARGS ####
#################################

datatypes <- c("RNAseq", "ATACseq", "H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3")

#############################################
## CREATE DIRECTORY STRUCTURE & LOAD DATA  ##
#############################################

## Set directory structure
folder = "1_DE_infection_modeling"

## Set directory 3 steps above script.
setwd('../../../')

for (i in 1:length(datatypes)){

  data <- datatypes[i]
  print(data)

  system(paste0("mkdir -p Outputs/",folder,"/", data, "_Mock_residuals"))
  system(paste0("mkdir -p Outputs/", folder,"/", data,"_Mock_voom"))

  ## Load metadata file
  meta_data = read.table(paste0("Inputs/metadata/", data, "_metadata.txt"))
  ## Remove IPSC samples which are not used in analysis.
  meta_data <- meta_data[!grepl("EU122", meta_data$Genotyping_ID),]
  #remove data type string from chipseq metadata
  rownames(meta_data)<- sub("_H3K27ac", "", rownames(meta_data))
  rownames(meta_data)<- sub("_H3K4me1", "", rownames(meta_data))
  rownames(meta_data)<- sub("_H3K27me3", "", rownames(meta_data))
  rownames(meta_data)<- sub("_H3K4me3", "", rownames(meta_data))

  if(data=="RNAseq"){
    reads = read.table(paste0("Inputs/counts_matrices/", data, "_filtered.counts.txt"), header=TRUE, row.names=1)
  } else {
    #Load read counts (row 1 as rownames)
    reads = read.table(paste0("Inputs/counts_matrices/", data, "_counts_mock.txt"), header=TRUE, row.names=1)
  }

  ## Set output directories
  voom_dir <- paste0("Outputs/", folder,"/", data,"_Mock_voom")
  results_dir <- paste0("Outputs/", folder,"/", data,"_Mock_results")
  residuals_dir <- paste0( "Outputs/", folder, "/", data, "_Mock_residuals")

  #################################
  ## PREPARE FILES FOR ANALYSIS  ##
  #################################

  ## Factorize certain columns from meta data
  meta_data$Batch <- as.factor(meta_data$Batch)
  meta_data$Genotyping_ID <- as.factor(meta_data$Genotyping_ID)

  ## Set original names in case reordering is necessary.
  reorder_names <- rownames(meta_data)
  length(reorder_names)

  ## Subset reads to include only the samples in metadata file.
  if(length(reorder_names) == dim(reads)[2]){
    reads <- reads[reorder_names]
  ## If the values do not match, remove those rows from the reads counts file (this will remove Mock samples)
  }else{
    correct_names <- colnames(reads)
    meta_data <- meta_data[rownames(meta_data) %in% correct_names,]
    reorder_names <- rownames(meta_data)
    reads <- reads[reorder_names]
  }

  ## Check to make sure reads and meta match 
  length(reorder_names) == dim(reads)[2]


  #################################
  ## VOOM NORMALIZE READ COUNTS  ##
  #################################

  ## DGEList simple takes reads and converts to a DGE object format which is needed in downstream functions.
  dge <- DGEList(counts = reads)
  ## CalcNormFactors does normalization by library size (default is TMM)
  dge <- calcNormFactors(dge)
  ## Set design matrix
  design = model.matrix(~0 + Batch, data = meta_data)
  ## Remove any columns that are all 0s
  design <- design[, colSums(design != 0) > 0]
  v <- voom(dge, design, plot = FALSE)
  ## Apply mean variance weights and design matrix to phenotype data.
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  #################################
  ##### GET & SAVE RESIDUALS ######
  #################################

  ## Get residuals to regress out batch effect
  residuals <- residuals.MArrayLM(object = fit, v)
  ## residuals should be the same dim as the expression matrix
  length(v$E) == length(residuals)

  ## Calculate average batch effect for each gene
  avg_batch_effect <- rowMeans(fit$coefficients)

  ## Add the average batch effect back into residuals
  corrected_expression <- apply(residuals,2,function(x){x + avg_batch_effect})
  weights <- v$weights
  colnames(weights) <- colnames(corrected_expression)
  rownames(weights) <- rownames(corrected_expression)

  ## Write batch-corrected expression and weights
  write.table(corrected_expression, paste0(residuals_dir,"/corrected_expression_MOCK.txt"), quote = FALSE, sep = ",")
  write.table(weights, paste0(residuals_dir,"/weights_MOCK.txt"), quote = FALSE, sep = ",")
}