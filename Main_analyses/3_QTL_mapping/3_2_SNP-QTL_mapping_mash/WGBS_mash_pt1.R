## part 1 of the WGBS mash pipeline extracts and saves the data in the proper format to run mash.
## part 2 of the pipeline applies mash to the data and saves results.

# Load libraries
library(ashr)
library(mashr)
library(gplots)
library(viridis)
library(ggplot2)
library(corrplot)
library(rmeta)
library(dplyr)
set.seed(321)

######################
## Set data to WGBS ##
######################

## Declare data type
data <- "WGBS"


##############################################
## Create directory structure and load data ##
##############################################

folder = "3_QTL_mapping"
## Set directory 3 steps above script.
setwd('../../../')

## Create directory structure to save outputs.
system(paste0("mkdir -p Outputs/",folder,"/SNP-QTL_mash/", data, "/"))
out_dir <- paste0("Outputs/",folder,"/SNP-QTL_mash/", data, "/")


## Load and create mash inputs
results_original_NI =read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/NI/", data, "_NI/raw_results/result_original.txt"),header=TRUE)
results_original_Flu =read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/Flu/", data, "_Flu/raw_results/result_original.txt"),header=TRUE)

## extract betas and se of betas from matrix eqtls for all snp-gene pairs
#betas is results_original$beta
#se of betas
results_original_NI$beta_se = results_original_NI$beta/results_original_NI$statistic
results_original_Flu$beta_se = results_original_Flu$beta/results_original_Flu$statistic


##################
## For ALL SNPs ##
##################

## Sort first by snps, then by gene.
FLU_beta_SE_matrix_ALL <- arrange(results_original_Flu, snps, gene)
NI_beta_SE_matrix_ALL <- arrange(results_original_NI, snps, gene)

## Will be true if both dfs are in the same order.
all(FLU_beta_SE_matrix_ALL$snps == NI_beta_SE_matrix_ALL$snps)
all(FLU_beta_SE_matrix_ALL$gene == NI_beta_SE_matrix_ALL$gene)

## Bind and rename column ids
all_conditions <- cbind(FLU_beta_SE_matrix_ALL, NI_beta_SE_matrix_ALL)
colnames(all_conditions) <- c("snps_flu", "gene_flu", "statistic_flu", "pvalue_flu", "FDR_flu", "flu_beta", "flu_beta_ses",
                               "snps_NI", "gene_NI", "statistic_NI", "pvalue_NI", "FDR_NI", "NI_beta", "NI_beta_ses")

# order by gene and by minimum pvalue between conditions.
all_conditions <- arrange(all_conditions, gene_flu, pmin(pvalue_flu, pvalue_NI))
## select best snp.
event.bestQTL<-all_conditions[!duplicated(all_conditions$gene_flu),]
## set rownames
rownames(event.bestQTL) <- event.bestQTL$gene_flu
## Replace any rows that have NAs for ses of betas. This means that these sites are invariable in one condition (will only be true for WGBS data).
all_conditions[is.na(all_conditions)] <- 0
dim(all_conditions)[1]


## Select betas and ses for ALL tests
betas_all <- dplyr::select(all_conditions, "flu_beta", "NI_beta")
betas_ses_all <- dplyr::select(all_conditions, "flu_beta_ses", "NI_beta_ses")

## Select betas and ses for BEST SNPs
betas_strong <- dplyr::select(event.bestQTL, "flu_beta", "NI_beta")
betas_ses_strong <- dplyr::select(event.bestQTL, "flu_beta_ses", "NI_beta_ses")


###################################################################
## GET RANDOM SUBSET OF ALL TESTS AND SET CORRELATION STRUCTURE  ##
###################################################################

## Get random subset of all tests performed
random.subset = sample(1:nrow(as.matrix(betas_all)), 200000)

## Estimate null correlation structure
data.temp = mash_set_data(as.matrix(as.matrix(betas_all)[random.subset,]), as.matrix(as.matrix(betas_ses_all)[random.subset,]), zero_Bhat_Shat_reset=2.22044604925031e-16)
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

## Transform inputs to mashr format
data_RANDOM = mash_set_data(as.matrix(betas_all[random.subset, ]), as.matrix(betas_ses_all[random.subset, ]), V = Vhat, zero_Bhat_Shat_reset=2.22044604925031e-16)
data_STRONG = mash_set_data(as.matrix(betas_strong), as.matrix(betas_ses_strong), V = Vhat, zero_Bhat_Shat_reset=2.22044604925031e-16)

## Set up the data-driven covariance matrix
U.pca = cov_pca(data_STRONG, 2)
## Extreme deconvolution
U.ed = cov_ed(data_STRONG, U.pca)
## Set up the canonical covariance matrix
U.c = cov_canonical(data_RANDOM)


###################################################
## Fit mash model (estimate mixture proportions) ##
###################################################

m = mash(data_RANDOM, Ulist = c(U.ed, U.c), outputlevel = 1)

## save this so you can load it in to the other parts
save(m, file=paste0(out_dir, "mash_random", ".Rdata"))


## save inputs for best matrix eqtl 
data_STRONG_p1 = mash_set_data(as.matrix(betas_strong[1:1000000, ]), as.matrix(betas_ses_strong[1:1000000, ]), V = Vhat, zero_Bhat_Shat_reset=2.22044604925031e-16)
data_STRONG_p2 = mash_set_data(as.matrix(betas_strong[1000001:2000000, ]), as.matrix(betas_ses_strong[1000001:2000000, ]), V = Vhat, zero_Bhat_Shat_reset=2.22044604925031e-16)
data_STRONG_p3 = mash_set_data(as.matrix(betas_strong[2000001:3000000, ]), as.matrix(betas_ses_strong[2000001:3000000, ]), V = Vhat, zero_Bhat_Shat_reset=2.22044604925031e-16)
data_STRONG_p4 = mash_set_data(as.matrix(betas_strong[3000001:4000000, ]), as.matrix(betas_ses_strong[3000001:4000000, ]), V = Vhat, zero_Bhat_Shat_reset=2.22044604925031e-16)
data_STRONG_p5 = mash_set_data(as.matrix(betas_strong[4000001:5000000, ]), as.matrix(betas_ses_strong[4000001:5000000, ]), V = Vhat, zero_Bhat_Shat_reset=2.22044604925031e-16)
data_STRONG_p6 = mash_set_data(as.matrix(betas_strong[5000001:6000000, ]), as.matrix(betas_ses_strong[5000001:6000000, ]), V = Vhat, zero_Bhat_Shat_reset=2.22044604925031e-16)
data_STRONG_p7 = mash_set_data(as.matrix(betas_strong[6000001:7000000, ]), as.matrix(betas_ses_strong[6000001:7000000, ]), V = Vhat, zero_Bhat_Shat_reset=2.22044604925031e-16)
data_STRONG_p8 = mash_set_data(as.matrix(betas_strong[7000001:nrow(betas_strong), ]), as.matrix(betas_ses_strong[7000001:nrow(betas_strong), ]), V = Vhat, zero_Bhat_Shat_reset=2.22044604925031e-16)

save(data_STRONG_p1, file=paste(out_dir, "data_STRONG_p1.Rdata"))
save(data_STRONG_p2, file=paste(out_dir, "data_STRONG_p2.Rdata"))
save(data_STRONG_p3, file=paste(out_dir, "data_STRONG_p3.Rdata"))
save(data_STRONG_p4, file=paste(out_dir, "data_STRONG_p4.Rdata"))
save(data_STRONG_p5, file=paste(out_dir, "data_STRONG_p5.Rdata"))
save(data_STRONG_p6, file=paste(out_dir, "data_STRONG_p6.Rdata"))
save(data_STRONG_p7, file=paste(out_dir, "data_STRONG_p7.Rdata"))
save(data_STRONG_p8, file=paste(out_dir, "data_STRONG_p8.Rdata"))

