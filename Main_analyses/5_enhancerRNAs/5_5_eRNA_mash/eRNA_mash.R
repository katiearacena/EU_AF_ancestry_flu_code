# Load libraries
library(ashr)
library(mashr)
library(gplots)
library(viridis)
library(ggplot2)
library(corrplot)
library(rmeta)
library(dplyr)
library(data.table)
set.seed(321)

############################
## Load command line args ##
############################

data <- "eRNA"

##############################################
## Create directory structure and load data ##
##############################################

folder = "5_enhancerRNAs"
## Set directory 3 steps above script.
setwd('../../../')

## Create directory structure to save outputs.
system(paste0("mkdir -p Outputs/",folder,"/SNP-QTL_mash/", data, "/"))
out_dir <- paste0("Outputs/",folder,"/SNP-QTL_mash/", data, "/")

## Load and create mash inputs
results_original_NI =as.data.frame(fread(paste0("Outputs/5_enhancerRNAs/SNP-QTL_mapping/",data,"/NI/3expPCs_NI/raw_results/result_original.txt")))
results_original_Flu =as.data.frame(fread(paste0("Outputs/5_enhancerRNAs/SNP-QTL_mapping/",data,"/Flu/4expPCs_Flu/raw_results/result_original.txt")))

results_original_NI$V1 <- NULL
results_original_Flu$V1 <- NULL
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

## use the fit from the previous run of mash by specifying g=get_fitted_g(m), fixg=TRUE to compute posterior summaries for any subset of tests
m2 = mash(data_STRONG, g = get_fitted_g(m), fixg = TRUE)


#################################
## Extract Posterior Summaries ##
#################################

# Get effects that are “significant”, which  means they have lfsr less than .10 in at least one condition
thresh <- .10
m.pairwise_PM <- get_pairwise_sharing(m2, lfsr_thresh=thresh, factor = 0.5)
m.sig <- get_significant_results(m2, thresh=thresh, sig_fn=get_lfdr)


################################################
## Identify condition-specific and shared QTL ##
################################################

#write this as an output and place in summary directory. ***
sink(paste0(out_dir, "summary_info.txt"))

print("how many are significant in at least one condition")
print(length(m.sig))
print("how many are significant in just Flu")
length(get_significant_results(m2, thresh=thresh, sig_fn=get_lfsr, conditions=1))
print("how many are significant in just NI")
length(get_significant_results(m2, thresh=thresh, sig_fn=get_lfsr, conditions=2))

#finish printing to file.
sink()

## write posterior outs
write.table(get_lfsr(m2), paste0(out_dir, "lfsr_output.txt"), quote = FALSE)
write.table(get_lfdr(m2), paste0(out_dir, "lfdr_output.txt"), quote = FALSE)
write.table(get_pm(m2), paste0(out_dir, "posteriorMeans.txt"), quote = FALSE)
write.table(get_psd(m2), paste0(out_dir, "posteriorStandardDevs.txt"), quote = FALSE)

## save R object
save(m2, file=paste0(out_dir, "mash_output.Rdata"))

## Write results tables with lfdr - local FDR
lfdr <- as.data.frame(get_lfdr(m2))
colnames(lfdr) <- c("lfdr_FLU", "lfdr_NI")

results_lfdr <- merge(event.bestQTL, lfdr, by=0)
rownames(results_lfdr) <- results_lfdr$Row.names
results_lfdr$Row.names <- NULL

## also merge lfsr - local false sign rate (more conservative).
lfsr <- as.data.frame(get_lfsr(m2))
colnames(lfsr) <- c("lfsr_FLU", "lfsr_NI")

results_lfdr_lfsr <- merge(results_lfdr, lfsr, by=0)
rownames(results_lfdr_lfsr) <- results_lfdr_lfsr$Row.names
results_lfdr_lfsr$Row.names <- NULL

# Write best snps with lfdr and lfsr columns added.
write.table(results_lfdr_lfsr, paste0(out_dir, "results_BOTH.mash.txt"))