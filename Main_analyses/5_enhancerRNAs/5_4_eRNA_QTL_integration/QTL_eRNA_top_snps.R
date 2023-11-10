### To integrate QTL across datatypes and collect ALL snps

## Load required libraries.
library(dplyr)
library(parallel)
library(data.table)
library(tidyverse)

args<-commandArgs(TRUE)

## Declare data type
QTL_type <- "eRNA"
condition <- args[1]

if(length(args)==0)
  {
    print("WARNING: No arguments supplied.")
  }
print(args)

#############################################################
## CREATE DIRECTORY STRUCTURE, SET FUNCTIONS, & LOAD DATA  ##
#############################################################

folder = "5_enhancerRNAs"
## Set directory 3 steps above script.
setwd('../../../')

## Create directory structure to save outputs.
system(paste0("mkdir -p Outputs/",folder,"/eRNA_SNP-QTL_integration_all_snps/", condition, "/"))
out_dir <- paste0("Outputs/",folder,"/eRNA_SNP-QTL_integration_all_snps/", condition, "/")

datatypes <- c("RNAseq", "ATACseq", "H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "WGBS")

best_snps<- list()
full_res <- list()

## load eRNA snps
if(condition=="NI"){
  eRNA_best <- as.data.frame(fread(paste0("Outputs/5_enhancerRNAs/SNP-QTL_mapping/eRNA/NI/3expPCs_NI/results_best_SNPs_with_qval.txt")))
  eRNA_full <- as.data.frame(fread(paste0("Outputs/5_enhancerRNAs/SNP-QTL_mapping/eRNA/NI/3expPCs_NI/raw_results/result_original.txt")))
} else if (condition=="Flu"){
  eRNA_best <- as.data.frame(fread(paste0("Outputs/5_enhancerRNAs/SNP-QTL_mapping/eRNA/Flu/4expPCs_Flu/results_best_SNPs_with_qval.txt")))
  eRNA_full <- as.data.frame(fread(paste0("Outputs/5_enhancerRNAs/SNP-QTL_mapping/eRNA/Flu/4expPCs_Flu/raw_results/result_original.txt")))
}

## load full results from QTL mapping
for (i in 1:length(datatypes)){
  data_i <- datatypes[i]
  print(data_i)
  full_res[[i]] <- as.data.frame(fread(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/", data_i, "/", condition, "/", data_i, "_", condition, "/raw_results/result_original.txt")))
  rownames(full_res[[i]]) <- full_res[[i]]$V1; full_res[[i]]$V1 <- NULL
  full_res[[i]] <- select(full_res[[i]], snps, gene, pvalue)
  colnames(full_res[[i]]) <- c("snps", paste0(data_i,"_feature"), paste0(data_i,"_pvalue"))
}

## get snps used in eRNA. remove those that are not to speed up code.
eRNA_snps_tested <- unique(eRNA_full$snps)

full_res_uniq_snps <- list()

## filter the full results to only keep unique snps. 
for (i in 1:length(datatypes)){
  data_i <- datatypes[i]
  print(data_i)
  res <- full_res[[i]]
  ## only keep snps that were tested in eRNA
  res <- res[res$snps %in% eRNA_snps_tested, ]
  ## sort based on snp, then pvalue
  res.sort<-res[order(res[,"snps"],res[,3]),]
  res.uniq.snps <-res.sort[!duplicated(res.sort$snps),]
  ## then only keep snp with strongest pvalue and record feature
  full_res_uniq_snps[[i]] <- res.uniq.snps
}


other_QTL_merged_df <- full_res_uniq_snps %>% reduce(full_join, by = "snps")

## also only select on top snp for eRNA
eRNA_full.sort<-eRNA_full[order(eRNA_full[,"snps"],eRNA_full[,3]),]
eRNA_full.uniq.snps <-eRNA_full.sort[!duplicated(eRNA_full.sort$snps),]
eRNA_full.uniq.snps$V1 <- NULL
## merge with eRNA dataframe
MERGED_df <- left_join(eRNA_full.uniq.snps, other_QTL_merged_df, by="snps")

## now that we have dataframe with snps as rows and pvalues for each datatype as columns
## now we need to determine if each of these is a significant QTL or not

## get qvalue threshold for eQTL. do this by thresholding by qvalue and then selecting the max pvalue that passes that threshold 
sig_eRNA_QTL_threshold <- as.vector(summary(eRNA_best[eRNA_best$qvalues < .10, ]$pvalue)[6])
thresholds <- readRDS("/project2/lbarreiro/users/katie/EU_AF_ancestry_flu/Make_figures/QTL_integration_thresholds.q.10.RData")
if(condition == "NI"){
  RNAseq_thresh <- thresholds[1]
  ATACseq_thresh <- thresholds[2]
  H3K27ac_thresh <- thresholds[3]
  H3K27me3_thresh <- thresholds[4]
  H3K4me1_thresh <- thresholds[5]
  H3K4me3_thresh <- thresholds[6]
  WGBS_thresh <-  thresholds[7]
} else if(condition == "Flu"){
  RNAseq_thresh <- thresholds[8]
  ATACseq_thresh <- thresholds[9]
  H3K27ac_thresh <- thresholds[10]
  H3K27me3_thresh <- thresholds[11]
  H3K4me1_thresh <- thresholds[12]
  H3K4me3_thresh <- thresholds[13]
  WGBS_thresh <-  thresholds[14]
}

MERGED_df$eRNAQTL_status <- ifelse(MERGED_df$pvalue < sig_eRNA_QTL_threshold, 1, 0)
MERGED_df$eQTL_status <- ifelse(MERGED_df$RNAseq_pvalue < RNAseq_thresh, 1, 0)
MERGED_df$caQTL_status <- ifelse(MERGED_df$ATACseq_pvalue < ATACseq_thresh, 1, 0)
MERGED_df$H3K27acQTL_status <- ifelse(MERGED_df$H3K27ac_pvalue < H3K27ac_thresh, 1, 0)
MERGED_df$H3K27meQTL_status <- ifelse(MERGED_df$H3K27me3_pvalue < H3K27me3_thresh, 1, 0)
MERGED_df$HK4me1QTL_status <- ifelse(MERGED_df$H3K4me1_pvalue < H3K4me1_thresh, 1, 0)
MERGED_df$H3K4me3QTL_status <- ifelse(MERGED_df$H3K4me3_pvalue < H3K4me3_thresh, 1, 0)
MERGED_df$meQTL_status <- ifelse(MERGED_df$WGBS_pvalue < WGBS_thresh, 1, 0)

print("saving output")
write.table(MERGED_df, paste0(out_dir, "merged_df_", condition, ".txt"))