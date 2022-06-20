########################################
### Load libraries and set functions ###
########################################

library(bsseq)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(qvalue)
library(DSS)
library(mashr)

data <- "WGBS"
# Decide whether to permute data (FALSE) or load previous permutations(TRUE)
reproduce <- FALSE
# Perform 10 permutations to use for FDR correction.
iterations <- 10

## Set directory structure
folder = "2_ancestry_effects_modeling"
## Set directory 3 steps above script.
setwd('../../../')
system(paste0("mkdir -p Outputs/",folder,"/PopDE_analysis/"))
out_dir <- paste0("Outputs/",folder,"/PopDE_analysis/")

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

##############################
### Load data and metadata ###
##############################

#Load metadata
meta_data = read.table(paste0("Inputs/metadata/", data, "_metadata.txt"))
## Load unsmoothed BSseq object
load("Inputs/counts_matrices/WGBS_filtered.counts.BSobj.RData")
BS.fit <- BSobj.fit.nolow

#subset to only include those samples for which there is a 1 in the "PopDEset" column 
meta_data=meta_data[which(meta_data$PopDE_set==1),]
#ensure that the Sample_ID are the rownames of the cols dataframe
rownames(meta_data)=meta_data$Sample
#ensure that the cols dataframe is ordered alphabetically (based on sample_id)
meta_data=meta_data[order(rownames(meta_data)),]
## Mean center age & refactorize 
meta_data=pretty_up_cols(meta_data)
## Set NI as first level so that it is used as the reference.
meta_data$Condition = factor(meta_data$Condition, levels=c("NI","Flu"))

## Check to make sure the number of samples in metadata and WGBS reads data matches.
dim(meta_data)[1] == dim(BS.fit)[2]
dim(BS.fit)[1]==7463164 # number of cpg sites there should be


##################################################
### Model DE within each condition (real data) ###
##################################################

## Set design matrix
design <- select(meta_data, Condition, Admixture, Age, Batch)
design

## run the linear model using DSS
lm = DMLfit.multiFactor(BS.fit, design=design, formula=~Condition + Condition:Admixture + Condition:Age + Condition:Batch)

## extract the admixture coefficients
NI_results = DMLtest.multiFactor(lm, coef="ConditionNI:Admixture")
FLU_results = DMLtest.multiFactor(lm, coef="ConditionFlu:Admixture")

## set rownames
rownames(NI_results) <- paste0(NI_results$chr, "_", NI_results$pos)
rownames(FLU_results) <- paste0(FLU_results$chr, "_", FLU_results$pos)

###########################################################
### Do permutation tests to define null model for PopDE ###
###########################################################

if(reproduce)
{
  ### to reproduce results previously generated load the permuted data.
  shuffled_pvals_adm_NI=read.table(paste0(out_dir, data, "/permuted_pvalues/popDE_NI.txt"))
  shuffled_pvals_adm_FLU=read.table(paste0(out_dir, data, "/permuted_pvalues/popDE_FLU.txt"))
  ##if reproduce is set to TRUE this part below will be skipped. If reproduce = FALSE then you the code below will permute the data.   (***need to review***)
}else{

  design_random=design

  for(iter in 1:iterations)
  {
    if(iter%%1==0)print(iter)

    design_random$Admixture=sample(design_random$Admixture)

    lm_rand = DMLfit.multiFactor(BS.fit, design=design_random, formula=~Condition + Condition:Admixture + Condition:Age + Condition:Batch)

    ## extract the admixture coefficients
    permuted_adm_NI = DMLtest.multiFactor(lm_rand, coef="ConditionNI:Admixture")
    permuted_adm_FLU = DMLtest.multiFactor(lm_rand, coef="ConditionFlu:Admixture")

    rownames(permuted_adm_NI) <- paste0(permuted_adm_NI$chr, "_", permuted_adm_NI$pos)
    rownames(permuted_adm_FLU) <- paste0(permuted_adm_FLU$chr, "_", permuted_adm_FLU$pos)

    if(iter==1)
    {
      shuffled_pvals_adm_FLU <-data.frame(x=permuted_adm_FLU$pvals)
      shuffled_pvals_adm_NI <-data.frame(x=permuted_adm_NI$pvals)

      rownames(shuffled_pvals_adm_FLU)=rownames(permuted_adm_FLU)
      rownames(shuffled_pvals_adm_NI)=rownames(permuted_adm_NI)

    } else {
      shuffled_pvals_adm_FLU <- cbind(shuffled_pvals_adm_FLU,x=permuted_adm_FLU$pvals)
      shuffled_pvals_adm_NI <- cbind(shuffled_pvals_adm_NI,x=permuted_adm_NI$pvals)

    }
  }
  system(paste0("mkdir -p ", out_dir, data, "/permuted_pvalues"))
  #Writing permutation tests p-values tables
  write.table(shuffled_pvals_adm_FLU,paste0(out_dir, data, "/permuted_pvalues/popDE_FLU.txt"))
  write.table(shuffled_pvals_adm_NI,paste0(out_dir, data, "/permuted_pvalues/popDE_NI.txt"))
}

# make results directories
system(paste0("mkdir -p ", out_dir, data, "/summary_stats/"))
system(paste0("mkdir -p ", out_dir, data, "/results/"))

##########################################################################
### Write results and add additional FDR correction methods to results ###
##########################################################################

## Remove NAs (will be same between real data and permutations)
FLU_results <- FLU_results[complete.cases(FLU_results), ]
NI_results <- NI_results[complete.cases(NI_results), ]

shuffled_pvals_adm_FLU <- shuffled_pvals_adm_FLU[complete.cases(shuffled_pvals_adm_FLU), ]
shuffled_pvals_adm_NI <- shuffled_pvals_adm_NI[complete.cases(shuffled_pvals_adm_NI), ]

dim(FLU_results)[1] == dim(shuffled_pvals_adm_FLU)[1]
dim(NI_results)[1] == dim(shuffled_pvals_adm_NI)[1]

##calculate BH adjusted p-values.
BH_Flu_pvals <- p.adjust(FLU_results$pvals, method= "BH", n=length(FLU_results$pvals))
BH_NI_pvals <- p.adjust(NI_results$pvals, method= "BH", n=length(NI_results$pvals))

##Bind BH adjusted p-values to results table.
results_FLU_BH=cbind(FLU_results, BH_Flu_pvals)
results_NI_BH=cbind(NI_results, BH_NI_pvals)

## Calculate qvalue
q_Flu_pvals <- empPvals(stat=-log10(FLU_results$pvals), stat0=-log10(as.matrix(shuffled_pvals_adm_FLU)), pool = T)
qvals_Flu <- qvalue(p=q_Flu_pvals)$qvalue
q_NI_pvals <- empPvals(stat=-log10(NI_results$pvals), stat0=-log10(as.matrix(shuffled_pvals_adm_NI)), pool = T)
qvals_NI <- qvalue(p=q_NI_pvals)$qvalue

## Bind qvalue to same tables.
results_FLU_BH_q=cbind(results_FLU_BH, qvals_Flu)
results_NI_BH_q=cbind(results_NI_BH, qvals_NI)

## Write results tables
write.table(results_FLU_BH_q, paste0(out_dir, data, "/results/popDE_FLU.txt"))
write.table(results_NI_BH_q, paste0(out_dir, data, "/results/popDE_NI.txt"))


##########################################################
### Identify condition-specific and shared PopDE genes ###
##########################################################

#write this as an output and place in summary directory.
sink(paste0(out_dir, data, "/summary_stats/summary_info.txt"))

#what are the PopDE genes in each condition using the different FDR thresholds?
print("Flu PopDE genes with BH <= 0.20")
print(summary(results_FLU_BH$BH_Flu_pvals <= 0.20))
print("NI PopDE genes with BH <= 0.20")
print(summary(results_NI_BH$BH_NI_pvals <= 0.20))
print("Flu PopDE genes with q <= 0.20")
print(summary(results_FLU_BH_q$qvals_Flu <= 0.20))
print("NI PopDE genes with q <= 0.20")
print(summary(results_NI_BH_q$qvals_NI <= 0.20))

#finish printing to file.
sink()


#####################
###### RUN MASH #####
#####################

#create the directory structure for the output files.
system(paste0("mkdir -p ", out_dir, data, "/mash_results"))
system(paste0("mkdir -p ", out_dir, data, "/mash_summary_stats"))

############################################
######### Step 1: Read in the data #########
############################################
## To run mash you need data consisting of a matrix of effects (Bhat) and a matrix of standard errors (Shat),
## for 𝐽 effects (rows) in 𝑅 conditions (columns).

#merge NI and FLU results (post qvalue)
results_FLU_BH_q$ID <- rownames(results_FLU_BH_q)
results_NI_BH_q$ID <- rownames(results_NI_BH_q)

both <- full_join(results_NI_BH_q, results_FLU_BH_q, by="ID")
rownames(both) <- both$ID; both$ID <- NULL

## extract Wald statistics 
betas_admixture <- select(both, "stat.x", "stat.y")
colnames(betas_admixture) <- c("ConditionNI:Admixture", "ConditionFlu:Admixture")

## Set matrix of standard errors to be all 1s.
SE_admixture <- matrix(1, dim(betas_admixture)[1], dim(betas_admixture)[2])
colnames(SE_admixture) <- colnames(betas_admixture)
rownames(SE_admixture) <- rownames(SE_admixture)

mash_data = mash_set_data(as.matrix(betas_admixture), as.matrix(SE_admixture))


############################################
## Step 2: set up the covariance matrices ##
############################################

## Get random subset of all tests performed
random.subset = sample(1:nrow(as.matrix(betas_admixture)), 200000)

## Estimate null correlation structure
data.temp = mash_set_data(as.matrix(as.matrix(betas_admixture)[random.subset,]), as.matrix(as.matrix(SE_admixture)[random.subset,]))
Vhat = estimate_null_correlation_simple(data.temp)
rm(data.temp)

data_RANDOM = mash_set_data(as.matrix(betas_admixture[random.subset, ]), as.matrix(SE_admixture[random.subset, ]), V = Vhat)
data_STRONG = mash_set_data(as.matrix(betas_admixture), as.matrix(SE_admixture), V = Vhat)

## use the  “canonical” covariance matries
U.c = cov_canonical(mash_data)
print(names(U.c))

## Set up the data-driven covariance matrix
U.pca = cov_pca(data_STRONG, 2)
## Extreme deconvolution
U.ed = cov_ed(data_STRONG, U.pca)
## Set up the canonical covariance matrix
U.c = cov_canonical(data_RANDOM)

############################################
########## Step 3: fit the model ###########
############################################

## run model with both canonical and data-driven cov matrices using random subset of tests.
## this is where mash learns that many tests are null and corrects for it.
m = mash(data_RANDOM, Ulist = c(U.ed, U.c), outputlevel = 1)

## use the fit from the previous run of mash by specifying g=get_fitted_g(m), fixg=TRUE to compute posterior summaries for any subset of tests
m2 = mash(data_STRONG, g = get_fitted_g(m), fixg = TRUE)

############################################
### Step 4: Extract Posterior Summaries ####
############################################
# Get effects that are “significant”, which in our definition means they have lfsr less than .10 in at least one condition
thresh <- .10
m.pairwise_PM <- get_pairwise_sharing(m2, lfsr_thresh=thresh, factor = 0.5)
m.sig <- get_significant_results(m2, thresh=thresh, sig_fn=get_lfdr)

############################################################
#### Identify condition-specific and shared PopDE genes ####
############################################################

#write this as an output and place in summary directory.
sink(paste0(out_dir, data, "/mash_summary_stats/summary_info.txt"))

print("how many are significant in at least one condition")
print(length(m.sig))
print("how many are significant in just Flu")
length(get_significant_results(m2, thresh=thresh, sig_fn=get_lfsr, conditions=1))
print("how many are significant in just NI")
length(get_significant_results(m2, thresh=thresh, sig_fn=get_lfsr, conditions=2))

#finish printing to file.
sink()

## write posterior outs
write.table(get_lfsr(m2), paste0(out_dir, data, "/mash_results/lfsr_output.txt"), quote = FALSE)
write.table(get_lfdr(m2), paste0(out_dir, data, "/mash_results/lfdr_output.txt"), quote = FALSE)
write.table(get_pm(m2), paste0(out_dir, data, "/mash_results/posteriorMeans.txt"), quote = FALSE)
write.table(get_psd(m2), paste0(out_dir, data, "/mash_results/posteriorStandardDevs.txt"), quote = FALSE)

## save R object
save(m2, file=paste0(out_dir, data, "/mash_results/mash_output.Rdata"))

## Write results tables with lfdr - local FDR
lfdr <- as.data.frame(get_lfdr(m2))
lfdr_FLU <- select(lfdr, "ConditionFlu:Admixture")
colnames(lfdr_FLU) <- c("lfdr")

results_FLU_BH_q$ID <- rownames(results_FLU_BH_q)
lfdr_FLU$ID <- rownames(lfdr_FLU)
results_FLU_lfdr <- full_join(results_FLU_BH_q, lfdr_FLU, by="ID")
rownames(results_FLU_lfdr) <- results_FLU_lfdr$ID

## also merge lfsr - local false sign rate (more conservative).
lfsr <- as.data.frame(get_lfsr(m2))
lfsr_FLU <- select(lfsr, "ConditionFlu:Admixture")
colnames(lfsr_FLU) <- c("lfsr")

lfsr_FLU$ID <- rownames(lfsr_FLU)
results_FLU_lfdr_lfsr <- full_join(results_FLU_lfdr, lfsr_FLU, by="ID")
rownames(results_FLU_lfdr_lfsr) <- results_FLU_lfdr_lfsr$ID; results_FLU_lfdr_lfsr$ID <- NULL


## now for NI
lfdr_NI<- select(lfdr, "ConditionNI:Admixture")
colnames(lfdr_NI) <- c("lfdr")

results_NI_BH_q$ID <- rownames(results_NI_BH_q)
lfdr_NI$ID <- rownames(lfdr_NI)
results_NI_lfdr <- full_join(results_NI_BH_q, lfdr_NI, by="ID")
rownames(results_NI_lfdr) <- results_NI_lfdr$ID

## also	merge lfsr - local false sign rate (more conservative).
lfsr_NI <- select(lfsr, "ConditionNI:Admixture")
colnames(lfsr_NI) <- c("lfsr")

lfsr_NI$ID <- rownames(lfsr_NI)
results_NI_lfdr_lfsr <- merge(results_NI_lfdr, lfsr_NI, by="ID")
rownames(results_NI_lfdr_lfsr) <- results_NI_lfdr_lfsr$ID; results_NI_lfdr_lfsr$ID <- NULL

# Write mash results
write.table(results_FLU_lfdr_lfsr, paste0(out_dir, data, "/mash_results/popDE_FLU.txt"))
write.table(results_NI_lfdr_lfsr, paste0(out_dir, data, "/mash_results/popDE_NI.txt"))



## Add average methylation counts and differences to results file
cts <- getMeth(BS.fit, type="raw")

AF_meta <- meta_data[meta_data$Ethnicity == "AF", ]
EU_meta <- meta_data[meta_data$Ethnicity == "EU", ]

AF_list <- AF_meta$Sample
EU_list <- EU_meta$Sample

NI_samples <- cts[,grepl("NI", colnames(cts))]
FLU_samples <- cts[,grepl("Flu", colnames(cts))]

AF_NI_samples <- NI_samples[,colnames(NI_samples) %in% AF_list]
AF_FLU_samples <- FLU_samples[,colnames(FLU_samples) %in% AF_list]

EU_NI_samples <- NI_samples[,colnames(NI_samples) %in% EU_list]
EU_FLU_samples <- FLU_samples[,colnames(FLU_samples) %in% EU_list]

##################################
## GET DIFF  B/W BINARIZED POPS ##
##################################

### *** NEED TO MAKE SURE WORKS WITH OUT PUT REVIEW### 
############
## FOR NI ## 
############

AF_NI_avg <- as.data.frame(rowMeans(AF_NI_samples, na.rm=TRUE))
colnames(AF_NI_avg) <- c("AF_mean")

EU_NI_avg <- as.data.frame(rowMeans(EU_NI_samples, na.rm=TRUE))
colnames(EU_NI_avg) <- c("EU_mean")

AF_NI_avg$feature_id <- rownames(AF_NI_avg)
EU_NI_avg$feature_id <- rownames(EU_NI_avg)

NI_group_means <- full_join(AF_NI_avg, EU_NI_avg, by="feature_id")
rownames(NI_group_means) <- NI_group_means$feature_id

results_NI_BH_q$feature_id <- rownames(results_NI_BH_q)
popDE_NI_with_means <-  left_join(results_NI_BH_q, NI_group_means, by="feature_id")
popDE_NI_with_means$diff <- popDE_NI_with_means$AF_mean - popDE_NI_with_means$EU_mean

## positive value == more meth in AF inds
## negative value == more meth in EU inds

#############
## FOR FLU ## 
#############

AF_FLU_avg <- as.data.frame(rowMeans(AF_FLU_samples, na.rm=TRUE))
colnames(AF_FLU_avg) <- c("AF_mean")

EU_FLU_avg <- as.data.frame(rowMeans(EU_FLU_samples, na.rm=TRUE))
colnames(EU_FLU_avg) <- c("EU_mean")

AF_FLU_avg$feature_id <- rownames(AF_FLU_avg)
EU_FLU_avg$feature_id <- rownames(EU_FLU_avg)

FLU_group_means <- full_join(AF_FLU_avg, EU_FLU_avg, by="feature_id")
rownames(FLU_group_means) <- FLU_group_means$feature_id

results_FLU_BH_q$feature_id <- rownames(results_FLU_BH_q)
popDE_FLU_with_means <-  left_join(results_FLU_BH_q, FLU_group_means, by="feature_id")
popDE_FLU_with_means$diff <- popDE_FLU_with_means$AF_mean - popDE_FLU_with_means$EU_mean

## positive value == more meth in AF inds
## negative value == more meth in EU inds

rownames(popDE_NI_with_means) <- popDE_NI_with_means$feature_id; popDE_NI_with_means$feature_id <- NULL
rownames(popDE_FLU_with_means) <- popDE_FLU_with_means$feature_id; popDE_FLU_with_means$feature_id <- NULL

## write
write.table(popDE_FLU_with_means, paste0(out_dir, data, "/results/popDE_FLU.txt"))
write.table(popDE_NI_with_means, paste0(out_dir, data, "/results/popDE_NI.txt"))