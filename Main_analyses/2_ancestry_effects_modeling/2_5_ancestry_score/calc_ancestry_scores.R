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
gene_list <- Reduce(union, list(gene_list_tfna, gene_list_inflam, gene_list_IL2, gene_list_IL6, gene_list_ifng))

#############################################
## CREATE DIRECTORY STRUCTURE & LOAD DATA  ##
#############################################

## Set directory structure
folder = "2_ancestry_effects_modeling"
## Set directory 3 steps above script.
setwd('../../../')

system(paste0("mkdir -p Outputs/",folder,"/ancestry_score/"))
out_dir <- paste0("Outputs/",folder,"/ancestry_score/")

datatypes <- c("RNAseq","ATACseq","H3K27ac","H3K27me3","H3K4me1", "H3K4me3", "WGBS")
results_with_genes <- list()
counts_matrices <- list()
res <- list()
meta_data <- list()


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

  ## Get gene-peak pairs info
  if (data_type_i=="RNAseq"){
    results_with_genes[[i]] <- read.table(paste0("Outputs/1_DE_infection_modeling/GSEA/",data_type_i,"/resultsALL_with_genes.txt"), header=TRUE, sep=" ")
  }else {
    results_with_genes[[i]] <- read.table(paste0("Outputs/1_DE_infection_modeling/GSEA/",data_type_i,"/resultsALL_with_genes.txt"), header=TRUE, sep=" ")
    rownames(results_with_genes[[i]]) <- results_with_genes[[i]]$peak_ids
  }

  ## Get age and batch corrected counts matrix
  if (data_type_i =="WGBS"){
    load("Inputs/counts_matrices/WGBS_filtered.counts.BSobj.RData")
    cts <- getMeth(BSobj.fit.nolow, type="raw")
    counts_matrices[[i]] <- data.frame(cts)

  } else {
    counts_matrices[[i]] <- read.table(paste0("Outputs/2_ancestry_effects_modeling/batch_corrected_cts_matrices/", data_type_i, "_batch.age.corrected.txt"), header=TRUE, sep="\t")
  }

  ## next need to connect gene_id with counts matrix as a column.
  ## add feature id to counts_matrices and results_with_genes from rownames
  counts_matrices[[i]]$feature_id <- rownames(counts_matrices[[i]])
  results_with_genes[[i]]$feature_id <- rownames(results_with_genes[[i]])
  ## merge
  res[[i]] <- full_join(counts_matrices[[i]], results_with_genes[[i]], by="feature_id")

}

#############################
####### FOR IFN alpha #######
#############################

## subset on genes and take the average
## all genes/peaks not just best!
gene_values <- list()
NI_Flu_df <- list()
gene_values_popDE_sig <- list()
sig_testing <- list()

for(i in 1:length(datatypes)){
  
  data_type_i <- datatypes[i]
  
  ## load popDE hits 
  NI_file_name <- paste0("popDE_NI.txt")
  FLU_file_name <- paste0("popDE_FLU.txt")
  qval_working_dir <- paste0("Outputs/2_ancestry_effects_modeling/PopDE_analysis/", data_type_i,"/results/")
  
  NI_q<- read.table(paste0(qval_working_dir, NI_file_name))
  FLU_q <- read.table(paste0(qval_working_dir, FLU_file_name))
  
  NI_q$feature <- rownames(NI_q)
  FLU_q$feature <- rownames(FLU_q)
  
  popDE_res <- full_join(NI_q, FLU_q, by="feature")
  ## which features are popDE in either of the conditions
  sig_popDE_res <- popDE_res[popDE_res$qvals_NI < .10 | popDE_res$qvals_Flu < .10, ]
  sig_popDE_res <- sig_popDE_res$feature
  length(sig_popDE_res)
  
  sig_popDE_NI <- popDE_res[popDE_res$qvals_NI < .10, ]
  sig_popDE_NI <- sig_popDE_NI$feature
  length(sig_popDE_NI)
  
  sig_popDE_Flu <- popDE_res[popDE_res$qvals_Flu < .10, ]
  sig_popDE_Flu <- sig_popDE_Flu$feature
  length(sig_popDE_Flu)
  
  ##get all genes in hallmark gene set
  gene_values[[i]] <- res[[i]][which((res[[i]]$symbol) %in% gene_list_ifna),]
  
  ## further filter to only include those that are popDE
  gene_values_popDE_sig[[i]] <- gene_values[[i]][which(gene_values[[i]]$feature_id %in% sig_popDE_res), ]
  
  # match meta data
  meta_data_i <- meta_data[[i]]
  counts <- gene_values_popDE_sig[[i]] 
  rownames(counts) <- gene_values_popDE_sig[[i]]$feature_id
  
  #rownames(counts) <- counts$symbol
  reorder_names <- rownames(meta_data_i)
  counts <- counts[reorder_names]
  
  length(which(colnames(counts)!=rownames(meta_data_i)))
  
  ## scale by condition
  NI_samples <- rownames(meta_data_i[meta_data_i$Condition == "NI",])
  FLU_samples <- rownames(meta_data_i[meta_data_i$Condition == "Flu",])
  
  NI_to_plot <- counts[,colnames(counts) %in% NI_samples]
  FLU_to_plot <- counts[,colnames(counts) %in% FLU_samples]
  
  NI_scaled <- t(scale(t(NI_to_plot)))
  NI_means <- colMeans(NI_scaled, na.rm = TRUE)
  
  FLU_scaled <- t(scale(t(FLU_to_plot)))
  FLU_means <- colMeans(FLU_scaled, na.rm = TRUE)
  
  meta_data_i$Sample <- sub("_H3K27ac", "", meta_data_i$Sample)
  meta_data_i$Sample<- sub("_H3K4me1", "",meta_data_i$Sample)
  meta_data_i$Sample<- sub("_H3K27me3", "", meta_data_i$Sample)
  meta_data_i$Sample<- sub("_H3K4me3", "",meta_data_i$Sample)
  
  meta_data_i_NI <- meta_data_i[meta_data_i$Condition=="NI", ]
  meta_data_i_Flu <- meta_data_i[meta_data_i$Condition=="Flu", ]

  ## get NI df for plotting
  NI_df <- as.data.frame(NI_means)
  colnames(NI_df) <- c("ind_means")
  NI_df$Ancestry <- ifelse(rownames(NI_df) == meta_data_i_NI$Sample, meta_data_i_NI$Ethnicity, "NA")
  NI_df$Ancestry <- replace(NI_df$Ancestry, NI_df$Ancestry==1, "AF")
  NI_df$Ancestry <- replace(NI_df$Ancestry, NI_df$Ancestry==2, "EU")
  NI_df$Condition <- "NI"
  
  ## get Flu df for plotting
  FLU_df <- as.data.frame(FLU_means)
  colnames(FLU_df) <- c("ind_means")
  FLU_df$Ancestry <- ifelse(rownames(FLU_df) == meta_data_i_Flu$Sample, meta_data_i_Flu$Ethnicity, "NA")
  FLU_df$Ancestry <- replace(FLU_df$Ancestry, FLU_df$Ancestry==1, "AF")
  FLU_df$Ancestry <- replace(FLU_df$Ancestry, FLU_df$Ancestry==2, "EU")
  FLU_df$Condition <- "Flu"
  
  NI_Flu_df[[i]] <- rbind(NI_df, FLU_df)
  
  NI_Flu_df[[i]]$data <-  data_type_i 
  
  ## significance testing
  AF_NI_df <- NI_df[NI_df$Ancestry=="AF", ]
  EU_NI_df <- NI_df[NI_df$Ancestry=="EU", ]
  
  AF_FLU_df <- FLU_df[FLU_df$Ancestry=="AF", ]
  EU_FLU_df <- FLU_df[FLU_df$Ancestry=="EU", ]
  
  
  NI_pval <- wilcox.test(AF_NI_df$ind_means, y = EU_FLU_df$ind_means, alternative = c("two.sided"))$p.value
  FLU_pval <-wilcox.test(AF_FLU_df$ind_means, y = EU_FLU_df$ind_means, alternative = c("two.sided"))$p.value
  
  pvals <- data.frame(NI_pval, FLU_pval)
  rownames(pvals) <- data_type_i
  sig_testing[[i]] <- pvals
  
}

IFNalpha <- do.call(rbind.data.frame, NI_Flu_df)
write.table(IFNalpha, file= paste0(out_dir, "IFNalpha.df.txt"))

IFNalpha_sig_testing_res <- do.call(rbind.data.frame, sig_testing)
write.table(IFNalpha_sig_testing_res, file=paste0(out_dir,"IFNalpha.pvals.txt"))

#################################################
####### FOR 5 other inflammatory pathways #######
#################################################

gene_values <- list()
NI_Flu_df <- list()
gene_values_popDE_sig <- list()
sig_testing <- list()

for(i in 1:length(datatypes)){
  
  data_type_i <- datatypes[i]
  
  ## load popDE hits 
  NI_file_name <- paste0("popDE_NI.txt")
  FLU_file_name <- paste0("popDE_FLU.txt")
  qval_working_dir <- paste0("Outputs/2_ancestry_effects_modeling/PopDE_analysis/", data_type_i,"/results/")
  
  NI_q<- read.table(paste0(qval_working_dir, NI_file_name))
  FLU_q <- read.table(paste0(qval_working_dir, FLU_file_name))
  
  NI_q$feature <- rownames(NI_q)
  FLU_q$feature <- rownames(FLU_q)
  
  popDE_res <- full_join(NI_q, FLU_q, by="feature")
  ## which features are popDE in either of the conditions
  sig_popDE_res <- popDE_res[popDE_res$qvals_NI < .10 | popDE_res$qvals_Flu < .10, ]
  sig_popDE_res <- sig_popDE_res$feature
  length(sig_popDE_res)
  
  sig_popDE_NI <- popDE_res[popDE_res$qvals_NI < .10, ]
  sig_popDE_NI <- sig_popDE_NI$feature
  length(sig_popDE_NI)
  
  sig_popDE_Flu <- popDE_res[popDE_res$qvals_Flu < .10, ]
  sig_popDE_Flu <- sig_popDE_Flu$feature
  length(sig_popDE_Flu)
  
  ##get all genes in hallmark gene set
  gene_values[[i]] <- res[[i]][which((res[[i]]$symbol) %in% gene_list),]
  
  ## further filter to only include those that are popDE
  gene_values_popDE_sig[[i]] <- gene_values[[i]][which(gene_values[[i]]$feature_id %in% sig_popDE_res), ]
  
  # match meta data
  meta_data_i <- meta_data[[i]]
  counts <- gene_values_popDE_sig[[i]] 
  rownames(counts) <- gene_values_popDE_sig[[i]]$feature_id
  
  #rownames(counts) <- counts$symbol
  reorder_names <- rownames(meta_data_i)
  counts <- counts[reorder_names]
  
  length(which(colnames(counts)!=rownames(meta_data_i)))
  
  ## scale by condition
  NI_samples <- rownames(meta_data_i[meta_data_i$Condition == "NI",])
  FLU_samples <- rownames(meta_data_i[meta_data_i$Condition == "Flu",])
  
  NI_to_plot <- counts[,colnames(counts) %in% NI_samples]
  FLU_to_plot <- counts[,colnames(counts) %in% FLU_samples]
  
  NI_scaled <- t(scale(t(NI_to_plot)))
  NI_means <- colMeans(NI_scaled, na.rm = TRUE)
  
  FLU_scaled <- t(scale(t(FLU_to_plot)))
  FLU_means <- colMeans(FLU_scaled, na.rm = TRUE)
  
  meta_data_i$Sample <- sub("_H3K27ac", "", meta_data_i$Sample)
  meta_data_i$Sample<- sub("_H3K4me1", "",meta_data_i$Sample)
  meta_data_i$Sample<- sub("_H3K27me3", "", meta_data_i$Sample)
  meta_data_i$Sample<- sub("_H3K4me3", "",meta_data_i$Sample)
  
  meta_data_i_NI <- meta_data_i[meta_data_i$Condition=="NI", ]
  meta_data_i_Flu <- meta_data_i[meta_data_i$Condition=="Flu", ]

  ## get NI df for plotting
  NI_df <- as.data.frame(NI_means)
  colnames(NI_df) <- c("ind_means")
  NI_df$Ancestry <- ifelse(rownames(NI_df) == meta_data_i_NI$Sample, meta_data_i_NI$Ethnicity, "NA")
  NI_df$Ancestry <- replace(NI_df$Ancestry, NI_df$Ancestry==1, "AF")
  NI_df$Ancestry <- replace(NI_df$Ancestry, NI_df$Ancestry==2, "EU")
  NI_df$Condition <- "NI"
  
  ## get Flu df for plotting
  FLU_df <- as.data.frame(FLU_means)
  colnames(FLU_df) <- c("ind_means")
  FLU_df$Ancestry <- ifelse(rownames(FLU_df) == meta_data_i_Flu$Sample, meta_data_i_Flu$Ethnicity, "NA")
  FLU_df$Ancestry <- replace(FLU_df$Ancestry, FLU_df$Ancestry==1, "AF")
  FLU_df$Ancestry <- replace(FLU_df$Ancestry, FLU_df$Ancestry==2, "EU")
  FLU_df$Condition <- "Flu"
  
  NI_Flu_df[[i]] <- rbind(NI_df, FLU_df)
  
  NI_Flu_df[[i]]$data <-  data_type_i 
  
  ## significance testing
  AF_NI_df <- NI_df[NI_df$Ancestry=="AF", ]
  EU_NI_df <- NI_df[NI_df$Ancestry=="EU", ]
  
  AF_FLU_df <- FLU_df[FLU_df$Ancestry=="AF", ]
  EU_FLU_df <- FLU_df[FLU_df$Ancestry=="EU", ]
  
  
  NI_pval <- wilcox.test(AF_NI_df$ind_means, y = EU_FLU_df$ind_means, alternative = c("two.sided"))$p.value
  FLU_pval <-wilcox.test(AF_FLU_df$ind_means, y = EU_FLU_df$ind_means, alternative = c("two.sided"))$p.value
  
  pvals <- data.frame(NI_pval, FLU_pval)
  rownames(pvals) <- data_type_i
  sig_testing[[i]] <- pvals
  
}


five.inflam_pathways <- do.call(rbind.data.frame, NI_Flu_df)
write.table(five.inflam_pathways, file=paste0(out_dir,"five.inflam_pathways.df.txt"))

five.inflam_pathways_sig_testing_res <- do.call(rbind.data.frame, sig_testing)
write.table(five.inflam_pathways_sig_testing_res, file=paste0(out_dir,"five.inflam_pathways.pvals.txt"))

