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

datatypes <- c("RNAseq","ATACseq","H3K27ac","H3K27me3","H3K4me1", "H3K4me3")
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
    counts_matrices[[i]] <- read.table(paste0("Outputs/4_cisregression_modeling/cisreg_ancestry_score/", data_type_i, "/corrected_matrix_noNAs.txt"), header=TRUE, sep=" ")
  }

  ## next need to connect gene_id with counts matrix as a column.
  ## add feature id to counts_matrices and results_with_genes from rownames
  counts_matrices[[i]]$feature_id <- rownames(counts_matrices[[i]])
  results_with_genes[[i]]$feature_id <- rownames(results_with_genes[[i]])
  ## merge and only keep those that are in counts matrix
  res[[i]] <- left_join(counts_matrices[[i]], results_with_genes[[i]], by="feature_id")

}


## subset on genes and take the average
## all genes/peaks not just best!
gene_values <- list()
NI_Flu_df <- list()
gene_values_popDE_sig <- list()
sig_testing <- list()
SET_df <- list()

for(i in 1:length(datatypes)){

  data_type_i <- datatypes[i]
  print(data_type_i)

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

  ##  filter to only include those that are popDE in either of the conditions
  gene_values_popDE_sig[[i]] <- res[[i]][which(res[[i]]$feature_id %in% sig_popDE_res), ]

  if (data_type_i=="RNAseq"){
    gene_values_popDE_sig[[i]] <- dplyr::select(gene_values_popDE_sig[[i]], !c(feature_id, betas, t.statistic, pvalues, fdrs, qvals))
  } else if (data_type_i =="WGBS"){
    gene_values_popDE_sig[[i]] <- dplyr::select(gene_values_popDE_sig[[i]], !c(feature_id, chr, pos, stat, pvals, fdrs, qvals, NI_mean, Flu_mean, diff, peaks, gene_ids))
  } else {
    gene_values_popDE_sig[[i]] <- dplyr::select(gene_values_popDE_sig[[i]], !c(feature_id, betas, t.statistic, pvalues, fdrs, qvals, gene_ids, peak_ids))
  }

  X <- as.data.table(gene_values_popDE_sig[[i]])
  ## take mean of peaks that share the same gene symbol
  X <- X[,lapply(.SD,mean),"symbol"]
  gene_values_popDE_sig[[i]] <- as.data.frame(X)
  ## remove any NAs (peaks not connected to genes)
  gene_values_popDE_sig[[i]] <- gene_values_popDE_sig[[i]][!is.na( gene_values_popDE_sig[[i]]$symbol),]
  rownames(gene_values_popDE_sig[[i]]) <- gene_values_popDE_sig[[i]]$symbol; gene_values_popDE_sig[[i]]$symbol <- NULL

  X2 <- gene_values_popDE_sig[[i]]

  # match meta data
  meta_data_i <- meta_data[[i]]

  #rownames(counts) <- counts$symbol
  reorder_names <- rownames(meta_data_i)
  X2 <- X2[reorder_names]

  length(which(colnames(X2)!=rownames(meta_data_i)))

  ## scale by condition
  NI_samples <- rownames(meta_data_i[meta_data_i$Condition == "NI",])
  FLU_samples <- rownames(meta_data_i[meta_data_i$Condition == "Flu",])

  NI_X2 <- X2[,colnames(X2) %in% NI_samples]
  FLU_X2 <- X2[,colnames(X2) %in% FLU_samples]

  meta_data_i$Sample <- sub("_H3K27ac", "", meta_data_i$Sample)
  meta_data_i$Sample<- sub("_H3K4me1", "",meta_data_i$Sample)
  meta_data_i$Sample<- sub("_H3K27me3", "", meta_data_i$Sample)
  meta_data_i$Sample<- sub("_H3K4me3", "",meta_data_i$Sample)

  meta_data_i_NI <- meta_data_i[meta_data_i$Condition=="NI", ]
  meta_data_i_Flu <- meta_data_i[meta_data_i$Condition=="Flu", ]

  length(which(is.na(rowVars(as.matrix(NI_X2), na.rm=F) == "NA")=="TRUE"))
  length(which(is.na(rowVars(as.matrix(FLU_X2), na.rm=F) == "NA")=="TRUE"))

  NI_gsva <- gsva(as.matrix(NI_X2), gs2, verbose=TRUE)
  Flu_gsva <- gsva(as.matrix(FLU_X2), gs2, verbose=FALSE)


  ## reset lists for each data type in case some results are not available for some pathways
  SET_pvals <- list()

  ### create plotting df for every row/ gene set
  for (j in 1:dim(NI_gsva)[1]){

    SET_NAME <- rownames(NI_gsva)[j]

    ## extract row you want to plot
    NI_gsva_set <- NI_gsva[j, ]
    NI_gsva_set <- as.data.frame(NI_gsva_set)
    ## get NI df for plotting
    Flu_gsva_set <- Flu_gsva[j, ]
    Flu_gsva_set <- as.data.frame(Flu_gsva_set)


    NI_df <- as.data.frame(NI_gsva_set)
    colnames(NI_df) <- c("ind_means")
    NI_df$Ancestry <- ifelse(rownames(NI_df) == meta_data_i_NI$Sample, meta_data_i_NI$Ethnicity, "NA")
    NI_df$Ancestry <- replace(NI_df$Ancestry, NI_df$Ancestry==1, "AF")
    NI_df$Ancestry <- replace(NI_df$Ancestry, NI_df$Ancestry==2, "EU")
    NI_df$Condition <- "NI"

    ## get Flu df for plotting
    FLU_df <- as.data.frame(Flu_gsva_set)
    colnames(FLU_df) <- c("ind_means")
    FLU_df$Ancestry <- ifelse(rownames(FLU_df) == meta_data_i_Flu$Sample, meta_data_i_Flu$Ethnicity, "NA")
    FLU_df$Ancestry <- replace(FLU_df$Ancestry, FLU_df$Ancestry==1, "AF")
    FLU_df$Ancestry <- replace(FLU_df$Ancestry, FLU_df$Ancestry==2, "EU")
    FLU_df$Condition <- "Flu"

    SET_df[[j]] <- rbind(NI_df, FLU_df)
    SET_df[[j]]$set <- SET_NAME

    ## significance testing
    AF_NI_df <- NI_df[NI_df$Ancestry=="AF", ]
    EU_NI_df <- NI_df[NI_df$Ancestry=="EU", ]

    AF_FLU_df <- FLU_df[FLU_df$Ancestry=="AF", ]
    EU_FLU_df <- FLU_df[FLU_df$Ancestry=="EU", ]


    NI_pval <- wilcox.test(AF_NI_df$ind_means, y = EU_FLU_df$ind_means, alternative = c("two.sided"))$p.value
    FLU_pval <-wilcox.test(AF_FLU_df$ind_means, y = EU_FLU_df$ind_means, alternative = c("two.sided"))$p.value

    SET_pvals[[j]] <- data.frame(NI_pval, FLU_pval)
    SET_pvals[[j]]$set <- SET_NAME

    }

  NI_Flu_df[[i]] <- do.call(rbind.data.frame, SET_df)
  NI_Flu_df[[i]]$data <-  data_type_i

  sig_testing[[i]] <- do.call(rbind.data.frame, SET_pvals)
  sig_testing[[i]]$data <-  data_type_i
}

GSVEA_SETS_FOR_PLOTTING <- do.call(rbind.data.frame, NI_Flu_df)
write.table(GSVEA_SETS_FOR_PLOTTING, file= paste0(out_dir, "GSVEA_SETS_FOR_PLOTTING.txt"))

GSVEA_sig_testing_res <- do.call(rbind.data.frame, sig_testing)
write.table(GSVEA_sig_testing_res, file=paste0(out_dir,"GSVEA_SETS_sig_testing_results.txt"))