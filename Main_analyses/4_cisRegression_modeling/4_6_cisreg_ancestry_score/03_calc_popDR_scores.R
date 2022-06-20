## Load libraries
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(data.table)
library(gridExtra)
library(mgsub)
library(metap)
library(qusage)
library(msigdbr)
library(biomaRt)
library(org.Hs.eg.db)
library(bsseq)
library(ggridges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(tidyverse)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(GSVA)
library(data.table)


## Get human pathways
m_df = msigdbr(species = "Homo sapiens")
## Check the available collections and sub-collections
m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
## Retrieve human genes for the hallmark collection gene sets
m_df = msigdbr(species = "Homo sapiens", category = "H")
human.path.list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

## select pathway you want gene list for
gene_list_tfna <- human.path.list$HALLMARK_TNFA_SIGNALING_VIA_NFKB
gene_list_inflam <- human.path.list$HALLMARK_INFLAMMATORY_RESPONSE
gene_list_IL2 <- human.path.list$HALLMARK_IL2_STAT5_SIGNALING
gene_list_IL6 <- human.path.list$HALLMARK_IL6_JAK_STAT3_SIGNALING
gene_list_ifna<- human.path.list$HALLMARK_INTERFERON_ALPHA_RESPONSE
gene_list_ifng <- human.path.list$HALLMARK_INTERFERON_GAMMA_RESPONSE

## get union of 5 inflammatory pathways (excluding ifna)
combo1 <- Reduce(union, list(gene_list_tfna, gene_list_inflam, gene_list_IL2, gene_list_IL6, gene_list_ifng))
combo2 <- Reduce(union, list(gene_list_tfna, gene_list_inflam, gene_list_IL2, gene_list_IL6))

## create list of all gene sets you want to test

gs2 <- list(gene_list_tfna, gene_list_inflam, gene_list_IL2, gene_list_IL6, gene_list_ifna, gene_list_ifng, combo1, combo2)
names(gs2) <- c("TFNA", "INFLAM", "IL2", "IL6", "IFNA", "IFNG", "TFNA_INFLA_IL2_IL6_IFNG", "TFNA_INFLA_IL2_IL6")

#############################################
## CREATE DIRECTORY STRUCTURE & LOAD DATA  ##
#############################################

## Set directory structure
folder = "4_cisregression_modeling"
## Set directory 3 steps above script.
setwd('../../../')

system(paste0("mkdir -p Outputs/",folder,"/cisreg_ancestry_score/"))
out_dir <- paste0("Outputs/",folder,"/cisreg_ancestry_score/")

## only include datatypes with sig popDR features (exclude H3K27me3 and WGBS)
datatypes <- c("RNAseq","ATACseq","H3K27ac", "H3K4me3")
results_with_genes <- list()
counts_matrices <- list()
res <- list()
meta_data <- list()
FC <- list()

## Create FC function
##Note: FC will be NA if one of the conditions has an NA.
get_FC=function(exp_stim,exp_ref,cols_stim,cols_ref,cond_stim,cond_ref){

  FC=exp_stim

  cont=0

  cols_all=rbind(cols_ref,cols_stim)
  exp_all=cbind(exp_ref,exp_stim)

  cols_all$Genotyping_ID=factor(cols_all$Genotyping_ID)
  cols_all$Condition=factor(cols_all$Condition)

  for(i in 1:length(levels(cols_all$Genotyping_ID)))
  {
    indexes=which(cols_all$Genotyping_ID==levels(cols_all$Genotyping_ID)[i])
    ref_index=indexes[which(cols_all$Condition[indexes]==cond_ref)]
    stim_index=indexes[which(cols_all$Condition[indexes]==cond_stim)]

    if(length(ref_index)==1 & length(stim_index)==1)
    {
      cont=cont+1
      FC[,cont]=exp_all[,stim_index]-exp_all[,ref_index]
      colnames(FC)[cont]=levels(cols_all$Genotyping_ID)[i]
    }
  }

  FC=FC[,1:cont]
  return(FC)
}

## Get counts matrix
for(i in 1:length(datatypes)){

  data_type_i <- datatypes[i]
  print(data_type_i)
 
  ## Load metadata
  meta_data[[i]] <- read.table(paste0("Inputs/metadata/", data_type_i, "_metadata.txt"))

  #subset to only include those samples for which there is a 1 in the "PopDEset" column
  meta_data[[i]]=meta_data[[i]][which(meta_data[[i]]$PopDE_set==1),]
  #ensure that the Sample_ID are the rownames of the cols dataframe
  rownames(meta_data[[i]])=meta_data[[i]]$Sample
  #ensure that the cols dataframe is ordered alphabetically (based on sample_id)
  meta_data[[i]]=meta_data[[i]][order(rownames(meta_data[[i]])),]
  #remove data type string from chipseq metadata
  rownames(meta_data[[i]])<- sub("_H3K27ac", "", rownames(meta_data[[i]]))
  rownames(meta_data[[i]])<- sub("_H3K4me1", "", rownames(meta_data[[i]]))
  rownames(meta_data[[i]])<- sub("_H3K27me3", "", rownames(meta_data[[i]]))
  rownames(meta_data[[i]])<- sub("_H3K4me3", "", rownames(meta_data[[i]]))

  ## Create FC matrices with cisreg matrices
  counts_matrices[[i]] <- read.table(paste0("Outputs/4_cisregression_modeling/cisreg_ancestry_score/", data_type_i, "/corrected_matrix_noNAs.txt"), header=TRUE, sep=" ")

  ### Get the condition-specific reads matrices
  ##subset the NI samples only)
  reads_NI=counts_matrices[[i]][,which(meta_data[[i]]$Condition=="NI")]
  ##subset the Flu samples only)
  reads_FLU=counts_matrices[[i]][,which(meta_data[[i]]$Condition=="Flu")]

  ### Get condition-specific metadata tables
  cols_NI=meta_data[[i]][which(meta_data[[i]]$Condition=="NI"),]
  cols_FLU=meta_data[[i]][which(meta_data[[i]]$Condition=="Flu"),]

  #this should be 0. It is a check to ensure that all the metadata and read counts matrix have all the same samples.
  length(which(rownames(cols_NI)!=colnames(reads_NI)))
  length(which(rownames(cols_FLU)!=colnames(reads_FLU)))

  FC[[i]]=get_FC(exp_stim=reads_FLU,exp_ref=reads_NI,cols_stim=cols_FLU,cols_ref=cols_NI,cond_stim="Flu",cond_ref="NI")

  ## Get gene-peak pairs info
  if (data_type_i=="RNAseq"){
    results_with_genes[[i]] <- read.table(paste0("Outputs/1_DE_infection_modeling/GSEA/",data_type_i,"/resultsALL_with_genes.txt"), header=TRUE, sep=" ")
  }else {
    results_with_genes[[i]] <- read.table(paste0("Outputs/1_DE_infection_modeling/GSEA/",data_type_i,"/resultsALL_with_genes.txt"), header=TRUE, sep=" ")
    rownames(results_with_genes[[i]]) <- results_with_genes[[i]]$peak_ids
  }


  ## next need to connect gene_id with FC matrix as a column.
  ## add feature id to counts_matrices and results_with_genes from rownames
  FC[[i]]$feature_id <- rownames(FC[[i]])
  results_with_genes[[i]]$feature_id <- rownames(results_with_genes[[i]])
  ## merge
  res[[i]] <- left_join(FC[[i]], results_with_genes[[i]], by="feature_id")
  rownames(res[[i]]) <- res[[i]]$feature_id

}


## subset on genes and take the average
## all genes/peaks not just best!
gene_values <- list()
NI_Flu_df <- list()
gene_values_popDR_sig <- list()
sig_testing <- list()
SET_df <- list()

for(i in 1:length(datatypes)){

  data_type_i <- datatypes[i]
  print(data_type_i)

  ## load popDR hits
  file_name <- paste0("results_FLU.txt")
  qval_working_dir <- paste0("Outputs/2_ancestry_effects_modeling/PopDR_analysis/", data_type_i,"/results/")

  popDR <- read.table(paste0(qval_working_dir, file_name))
  popDR$feature <- rownames(popDR)

  ## which features are popDR
  sig_popDR_res <- popDR[popDR$qvals_FLU < .20, ]
  sig_popDR_res <- sig_popDR_res$feature
  length(sig_popDR_res)

  ##  filter to only include those that are popDR
  gene_values_popDR_sig[[i]] <- res[[i]][which(res[[i]]$feature_id %in% sig_popDR_res), ]

  if (data_type_i=="RNAseq"){
    gene_values_popDR_sig[[i]] <- dplyr::select(gene_values_popDR_sig[[i]], !c(feature_id, betas, t.statistic, pvalues, fdrs, qvals))
  } else {
    gene_values_popDR_sig[[i]] <- dplyr::select(gene_values_popDR_sig[[i]], !c(feature_id, betas, t.statistic, pvalues, fdrs, qvals, gene_ids, peak_ids))
  }

  X <- as.data.table(gene_values_popDR_sig[[i]])
  ## take mean of peaks that share the same gene symbol
  X <- X[,lapply(.SD,mean),"symbol"]
  gene_values_popDR_sig[[i]] <- as.data.frame(X)
  ## remove any NAs (peaks not connected to genes)
  gene_values_popDR_sig[[i]] <- gene_values_popDR_sig[[i]][!is.na( gene_values_popDR_sig[[i]]$symbol),]
  rownames(gene_values_popDR_sig[[i]]) <- gene_values_popDR_sig[[i]]$symbol; gene_values_popDR_sig[[i]]$symbol <- NULL

  X2 <- gene_values_popDR_sig[[i]]

  # match meta data
  meta_data_i <- meta_data[[i]]

  # get sample IDs from metadata
  meta_data_i_samples <- meta_data_i[meta_data_i$Condition=="NI", ]
  rownames(meta_data_i_samples) <- meta_data_i_samples$Genotyping_ID

  reorder_names <- rownames(meta_data_i_samples)

  X2 <- X2[reorder_names]

  length(which(colnames(X2)!=rownames(meta_data_i_samples)))

  ## scale
  samples <- rownames(meta_data_i_samples)

  FC_X2 <- X2[,colnames(X2) %in% samples]

  meta_data_i_samples$Sample <- sub("_H3K27ac", "", meta_data_i_samples$Sample)
  meta_data_i_samples$Sample<- sub("_H3K4me1", "",meta_data_i_samples$Sample)
  meta_data_i_samples$Sample<- sub("_H3K4me3", "",meta_data_i_samples$Sample)

  gsva <- gsva(as.matrix(FC_X2), gs2, verbose=FALSE)

  ## reset lists for each data type in case some results are not available for some pathways
  SET_pvals <- list()

  ### create plotting df for every row/ gene set
  for (j in 1:dim(gsva)[1]){

    SET_NAME <- rownames(gsva)[j]

    ## extract row you want to plot
    FC_gsva_set <- gsva[j, ]
    FC_gsva_set <- as.data.frame(FC_gsva_set)

    FC_df <- as.data.frame(FC_gsva_set)
    colnames(FC_df) <- c("ind_means")
    FC_df$Ancestry <- ifelse(rownames(FC_df) == meta_data_i_samples$Genotyping_ID, meta_data_i_samples$Ethnicity, "NA")
    FC_df$Ancestry <- replace(FC_df$Ancestry, FC_df$Ancestry==1, "AF")
    FC_df$Ancestry <- replace(FC_df$Ancestry, FC_df$Ancestry==2, "EU")

    SET_df[[j]] <- FC_df
    SET_df[[j]]$set <- SET_NAME

    ## significance testing
    AF_df <- FC_df[FC_df$Ancestry=="AF", ]
    EU_df <- FC_df[FC_df$Ancestry=="EU", ]

    FC_pval <- wilcox.test(AF_df$ind_means, y = EU_df$ind_means, alternative = c("two.sided"))$p.value

    SET_pvals[[j]] <- data.frame(FC_pval)
    SET_pvals[[j]]$set <- SET_NAME

    }

  NI_Flu_df[[i]] <- do.call(rbind.data.frame, SET_df)
  NI_Flu_df[[i]]$data <-  data_type_i

  sig_testing[[i]] <- do.call(rbind.data.frame, SET_pvals)
  sig_testing[[i]]$data <-  data_type_i
}

GSVEA_SETS_FOR_PLOTTING <- do.call(rbind.data.frame, NI_Flu_df)
write.table(GSVEA_SETS_FOR_PLOTTING, file= paste0(out_dir, "popDR_GSVEA_SETS_FOR_PLOTTING.txt"))

GSVEA_sig_testing_res <- do.call(rbind.data.frame, sig_testing)
write.table(GSVEA_sig_testing_res, file=paste0(out_dir,"popDR_GSVEA_SETS_sig_testing_results.txt"))