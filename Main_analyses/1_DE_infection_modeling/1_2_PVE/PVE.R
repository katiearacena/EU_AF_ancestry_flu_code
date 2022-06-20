## May have to subset and parallelize WGBS data in this analysis. 

# Load libraries
library(limma)
library(edgeR)
library(relaimpo)
library(robustbase)
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringi)
library(stringr)
library(ggpubr)
library(bsseq)

## Set directory structure
folder = "1_DE_infection_modeling"
## Set directory 2 steps above script.
setwd('../../../')
## set outputs directory
system(paste0("mkdir -p Outputs/",folder,"/PVE_analysis_results"))
OUTPUTS_dir <- paste0("Outputs/",folder,"/PVE_analysis_results")

## List data
datatypes <- c("RNAseq","ATACseq","H3K27ac","H3K27me3","H3K4me1", "H3K4me3", "WGBS")
rela_impo_outs <- data.frame(matrix(nrow = 2, ncol = length(datatypes)))

for(j in 1:length(datatypes)){

  data_type_i <- datatypes[j]
  print(paste0(data_type_i, " started"))

  ## Load metadata file
  meta_data = read.table(paste0("Inputs/metadata/", data_type_i, "_metadata.txt"))

  ## Factorize
  meta_data$Batch <- as.factor(meta_data$Batch)
  meta_data$Genotyping_ID <- as.factor(meta_data$Genotyping_ID)
  meta_data$Condition = factor(meta_data$Condition, levels=c("NI","Flu"))  
  #remove data type string from chipseq metadata
  rownames(meta_data)<- sub("_H3K27ac", "", rownames(meta_data))
  rownames(meta_data)<- sub("_H3K4me1", "", rownames(meta_data))
  rownames(meta_data)<- sub("_H3K27me3", "", rownames(meta_data))
  rownames(meta_data)<- sub("_H3K4me3", "", rownames(meta_data))

  ## Remove Mock and IPSC samples which are not used in analysis.
  if(data_type_i=="RNAseq"){
    meta_data <- meta_data[which(!meta_data$Condition=="Mock"),]
    meta_data <- meta_data[!grepl("EU122", meta_data$Genotyping_ID),]
  }else{
    meta_data <- meta_data[which(!meta_data$Condition=="Mock"),]
    meta_data <- meta_data[which(meta_data$PopDE_set=="1"),]
    meta_data <- meta_data[!grepl("EU122", meta_data$Genotyping_ID),]
    meta_data <- meta_data[!grepl("IPSCs", meta_data$Genotyping_ID),]
  }

  dim(meta_data)

  reorder_names <- rownames(meta_data)

  ## Read in methylation data 
  if(data_type_i=="WGBS"){
    load("Inputs/counts_matrices/WGBS_filtered.counts.BSobj.RData")

    BS.cov <- getCoverage(BSobj.fit.nolow)

    ## Filter based on coverage
    #To minimize noise in methylation estimates due to low-coverage data, we restricted analysis to CpG sites 
    # with coverage of ≥4 sequence reads in at least half of the samples in each condition. DSS accounts for low coverage
    # in the model so it is not necessary before. 

   keepLoci.ex <- which(rowSums(BS.cov[, BSobj.fit.nolow$Condition == "Flu"] >= 4) > 17 &
                       rowSums(BS.cov[, BSobj.fit.nolow$Condition == "NI"] >= 4) > 17 )
    length(keepLoci.ex)
    BSobj.fit.nolow <- BSobj.fit.nolow[keepLoci.ex,]
    BSobj.fit.nolow
    ## get cts matrix from bsseq object
    exp <- getMeth(BSobj.fit.nolow, type="raw")


    NI_samples <- exp[,grepl("NI", colnames(exp))]
    FLU_samples <- exp[,grepl("Flu", colnames(exp))]
       
    all.NA_NI <-  rownames(NI_samples[rowVars(NI_samples, na.rm=TRUE) == 0, ])
    all.NA_FLU <- rownames(FLU_samples[rowVars(FLU_samples, na.rm=TRUE) == 0, ])
    no.var <- rownames(exp[rowVars(exp, na.rm=TRUE) == 0, ])

    to_remove <- c(all.NA_NI, all.NA_FLU, no.var)
    to_remove_uniq <- unique(to_remove)
    length(to_remove_uniq)
    exp <- subset(exp, !(rownames(exp) %in% to_remove_uniq))

    ## should be == 0 if no mismatch
    length(which(colnames(exp)!=rownames(meta_data)))
    dim(exp)[1]

    ## model infection differential expression 
    exp = as.matrix(exp)
    gene_model = lm(exp[1,] ~ Genotyping_ID + Condition, data = meta_data)
    rel_impo = suppressWarnings(calc.relimp(gene_model))

    ## matrix to store
    importances = data.frame(matrix(nrow = nrow(exp), ncol = length(rel_impo$lmg)))
    rownames(importances) <- rownames(exp)
    colnames(importances) <- names(rel_impo$lmg)

    ## loop over all genes
    ## skip and place NA for any genes for which you cannot run relaimpo 
    ## (those with 0 variance in either condition or too many missing data points)
    for(i in 1:nrow(exp)){
        if(i%%1000 == 0) print(i)

        gene_model = lm(exp[i,] ~ Genotyping_ID + Condition, data = meta_data)
    
        out <- tryCatch(
            {
            rel_impo = suppressWarnings(calc.relimp(gene_model, type = "lmg"))
            importances[i,] = rel_impo$lmg
            }, 
            error=function(cond) {
            message(cond)
            message(" - SKIPPED")
            importances[i,] = "NA"
            }
            )  
        }

  ## if any other data type
    }else{
  ## read in corrected expression
    exp <- read.table(paste0("Outputs/", folder,"/", data_type_i, "_residuals/corrected_expression.txt"), header = TRUE, sep = ",", row.names = 1)
  ## read in weights
    weights <-  read.table(paste0("Outputs/", folder,"/", data_type_i, "_residuals/weights.txt"), header = TRUE, sep = ",", row.names = 1)

    ## Make sure order and dimensions match
    if(length(reorder_names) == dim(exp)[2]){
      exp <- exp[reorder_names]
      weights <- weights[reorder_names]
    }else{
      correct_names <- colnames(exp)
      meta_data <- meta_data[rownames(meta_data) %in% correct_names,]
      weights <- weights[,colnames(weights) %in% correct_names]
      reorder_names <- rownames(meta_data)
      exp <- exp[reorder_names]
      weights <- weights[reorder_names]
    }

    ## should all == 0 if no mismatch
    length(which(colnames(exp)!=rownames(meta_data)))
    length(which(colnames(weights)!=rownames(meta_data)))
    length(which(colnames(exp)!=colnames(weights)))

    ## model infection differential expression
    exp = as.matrix(exp)
    weights = as.matrix(weights)
    gene_model = lm(exp[1,] ~ Genotyping_ID + Condition, data = meta_data, weights = weights[1,])
    rel_impo = suppressWarnings(calc.relimp(gene_model))

    ## matrix to store
    importances = data.frame(matrix(nrow = nrow(exp), ncol = length(rel_impo$lmg)))
    rownames(importances) <- rownames(exp)
    colnames(importances) <- names(rel_impo$lmg)

    ## loop over all genes
    for(i in 1:nrow(exp))
    {
      if(i%%1000 == 0) print(i)
      gene_model = lm(exp[i,] ~ Genotyping_ID + Condition, data = meta_data, weights = weights[i,])
      rel_impo = suppressWarnings(calc.relimp(gene_model, type = "lmg"))
      importances[i,] = rel_impo$lmg
    }
  }

  write.table(importances, paste0(OUTPUTS_dir,"/importances_",data_type_i,".txt"), quote = FALSE)

  ## calculate col medians (median of all genes for each covariate -- collapses matrix)
  medians <- as.data.frame(colMedians(as.matrix(importances)))
  colnames(medians) <- c(paste0("medians_",data_type_i))

  rownames(rela_impo_outs) <- colnames(importances)
  colnames(rela_impo_outs)[j] <- colnames(medians)[1]
  rela_impo_outs[,j] <- medians[,1]

  print(paste0(data_type_i, " completed"))
}

write.table(rela_impo_outs, paste0(OUTPUTS_dir,"/rela_impo_outs.txt"), sep = ",", quote = FALSE)
