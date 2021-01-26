#load libraries
library(preprocessCore)
library(ggplot2)
library(limma)
library(edgeR)
library(statmod)
library(sva)
library(reshape2)
library(cobs)
library(dplyr)
library(qvalue)
library(corrplot)
library(ashr)
library(mashr)
library(rmeta)
library(grid)
library(gridExtra)

#################################
##### LOAD COMMAND LINE ARGS ####
#################################

## Load command line arguments
args<-commandArgs(TRUE)
## Datatype
data <- args[1]
## Reproduce permutations (TRUE or FALSE)
reproduce <- args[2]

if(length(args)==0)
{
  print("WARNING: No arguments supplied.")
}
print(args)


#############################################
## CREATE DIRECTORY STRUCTURE & LOAD DATA  ##
#############################################

## Set directory structure
folder = "2_ancestry_effects_modeling"
## Set directory 3 steps above script.
setwd('../../../')
system(paste0("mkdir -p Outputs/",folder,"/PopDE_analysis/"))
out_dir <- paste0("Outputs/",folder,"/PopDE_analysis/")

#Set inputs

#Load metadata
meta_data = read.table(paste0("Inputs/metadata/", data, "_metadata.txt"))
#Load read counts (row 1 as rownames)
reads = read.table(paste0("Outputs/", folder,"/batch_corrected_cts_matrices/",data, "_batch.age.corrected.txt"), header=TRUE, row.names=1)


########################################
### Format reads matrices & metadata ###
########################################

#subset to only include those samples for which there is a 1 in the "PopDEset" column
meta_data=meta_data[which(meta_data$PopDE_set==1),]
#ensure that the Sample_ID are the rownames of the cols dataframe
rownames(meta_data)=meta_data$Sample
#ensure that the cols dataframe is ordered alphabetically (based on sample_id)
meta_data=meta_data[order(rownames(meta_data)),]
#remove data type string from chipseq metadata
rownames(meta_data)<- sub("_H3K27ac", "", rownames(meta_data))
rownames(meta_data)<- sub("_H3K4me1", "", rownames(meta_data))
rownames(meta_data)<- sub("_H3K27me3", "", rownames(meta_data))
rownames(meta_data)<- sub("_H3K4me3", "", rownames(meta_data))

#make sure that all the columns of the meta data are also in the read counts matrix.
reads=reads[,which(colnames(reads) %in% rownames(meta_data))]
#order the reads count matrix.
reads=reads[,order(colnames(reads))]

#this should be 0. It is a check to ensure that all the metadata and read counts matrix have all the same samples.
length(which(rownames(meta_data)!=colnames(reads)))


#####################################################
#### Model DE within each condition (real data) #####
#####################################################

DE_WC=function(reads,cols){
  design = model.matrix(~0+ Condition + Condition:Admixture, data = cols)
  # Remove mock columns
  design <- design[, colSums(design != 0) > 0]
  fit <-lmFit(reads,design)
  fit <- eBayes(fit)
  return(list(fit,reads))
}

DE_BOTH=DE_WC(reads=reads,cols=meta_data)
fit_BOTH=DE_BOTH[[1]]
v_BOTH=DE_BOTH[[2]]


##############################################################
#### Do permutation tests to define null model for PopDE #####
##############################################################

# Perform 1000 permutations
iterations=1000

if(reproduce)
{
  ### to reproduce results previously generated load the permuted data.
  shuffled_pvals_adm_NI=read.table(paste0(out_dir, data,"/permuted_pvalues/popDE_NI.txt"))
  shuffled_pvals_adm_FLU=read.table(paste0(out_dir, data,"/permuted_pvalues/popDE_FLU.txt"))
  ##if reproduce is set to TRUE this part below will be skipped. If reproduce = FALSE then you the code below will permute the data.

}else{

  meta_data_random=meta_data

  for(iter in 1:iterations)
  {
    if(iter%%100==0)print(iter)
    ## permute admixture
    meta_data_random$Admixture=sample(meta_data_random$Admixture)

    DE_rand=DE_WC(reads=reads,cols=meta_data_random)

    permuted_adm_FLU=topTable(DE_rand[[1]], coef="ConditionFlu:Admixture", sort.by="none", adjust="BH",n=nrow(DE_rand[[1]]$coefficients))
    permuted_adm_NI=topTable(DE_rand[[1]], coef="ConditionNI:Admixture", sort.by="none", adjust="BH",n=nrow(DE_rand[[1]]$coefficients))

    if(iter==1)
    {
      shuffled_pvals_adm_FLU <-data.frame(x=permuted_adm_FLU$P.Value)
      shuffled_pvals_adm_NI <-data.frame(x=permuted_adm_NI$P.Value)

      rownames(shuffled_pvals_adm_FLU)=rownames(permuted_adm_FLU)
      rownames(shuffled_pvals_adm_NI)=rownames(permuted_adm_NI)

    } else {
      shuffled_pvals_adm_FLU <- cbind(shuffled_pvals_adm_FLU,x=permuted_adm_FLU$P.Value)
      shuffled_pvals_adm_NI <- cbind(shuffled_pvals_adm_NI,x=permuted_adm_NI$P.Value)

    }
  }
  system(paste0("mkdir -p ", out_dir, data, "/permuted_pvalues"))
  #Writing permutation tests p-values tables
  write.table(shuffled_pvals_adm_FLU,paste0(out_dir ,data,"/permuted_pvalues/popDE_FLU.txt"))
  write.table(shuffled_pvals_adm_NI,paste0(out_dir, data,"/permuted_pvalues/popDE_NI.txt"))

}

#create the directory structure for the output files.
system(paste0("mkdir -p ", out_dir, data, "/results"))
system(paste0("mkdir -p ", out_dir, data, "/summary_stats/"))

#create a data.frame of the results (and remove automatic p-value adjustment calculation)
results_FLU=topTable(fit_BOTH,coef="ConditionFlu:Admixture",number=nrow(fit_BOTH$coefficients))[,1:4]
results_NI=topTable(fit_BOTH,coef="ConditionNI:Admixture",number=nrow(fit_BOTH$coefficients))[,1:4]

##################################################
#### Write results and perform FDR correction ####
##################################################

##calculate BH adjusted p-values.
BH_Flu_pvals <- p.adjust(results_FLU$P.Value, method= "BH", n=length(results_FLU$P.Value))
BH_NI_pvals <- p.adjust(results_NI$P.Value, method= "BH", n=length(results_NI$P.Value))

##Bind BH adjusted p-values to results table.
results_FLU_BH=cbind(results_FLU, BH_Flu_pvals)
results_NI_BH=cbind(results_NI, BH_NI_pvals)

## Calculate qvalue
q_Flu_pvals <- empPvals(stat=-log10(results_FLU$P.Value), stat0=-log10(as.matrix(shuffled_pvals_adm_FLU)), pool = T)
qvals_Flu <- qvalue(p=q_Flu_pvals)$qvalue
q_NI_pvals <- empPvals(stat=-log10(results_NI$P.Value), stat0=-log10(as.matrix(shuffled_pvals_adm_NI)), pool = T)
qvals_NI <- qvalue(p=q_NI_pvals)$qvalue

## Bind qvalue to same tables.
results_FLU_BH_q=cbind(results_FLU_BH, qvals_Flu)
results_NI_BH_q=cbind(results_NI_BH, qvals_NI)

## Write results tables
write.table(results_FLU_BH_q, paste0(out_dir, data,"/results/popDE_FLU.txt"))
write.table(results_NI_BH_q, paste0(out_dir, data,"/results/popDE_NI.txt"))
## Write voomed read counts for plotting
write.table(v_BOTH$E, paste0(out_dir, data,"/results/voomed_counts_for_plotting.txt"))

###################################
#### Identify sig. PopDE genes ####
###################################

#write this as an output and place in summary directory.
sink(paste0(out_dir, data, "/summary_stats/summary_info.txt"))

#what are the PopDE genes in each condition using the different FDR thresholds?
print("Flu PopDE genes with BH <= 0.10")
print(summary(results_FLU_BH$BH_Flu_pvals <= 0.10))
print("NI PopDE genes with BH <= 0.10")
print(summary(results_NI_BH$BH_NI_pvals <= 0.10))
print("Flu PopDE genes with q <= 0.10")
print(summary(results_FLU_BH_q$qvals_Flu <= 0.10))
print("NI PopDE genes with q <= 0.10")
print(summary(results_NI_BH_q$qvals_NI <= 0.10))

#finish printing to file.
sink()

#####################
###### RUN MASH #####
#####################

#create the directory structure for the output files.
system(paste0("mkdir -p ", out_dir, data, "/mash_results"))
system(paste0("mkdir -p ", out_dir, data, "/mash_summary_stats"))

results_FLU=topTable(fit_BOTH,coef="ConditionFlu:Admixture",number=nrow(fit_BOTH$coefficients))
results_NI=topTable(fit_BOTH,coef="ConditionNI:Admixture",number=nrow(fit_BOTH$coefficients))

############################################
######### Step 1: Read in the data #########
############################################
## To run mash you need data consisting of a matrix of effects (Bhat) and a matrix of standard errors (Shat),
## for 𝐽 effects (rows) in 𝑅 conditions (columns).

betas <- as.data.frame(fit_BOTH$coefficients)
betas_admixture <- select(betas, "ConditionFlu:Admixture", "ConditionNI:Admixture")

## Get standard errors of betas 
SE <- as.data.frame(sqrt(fit_BOTH$s2.post) * fit_BOTH$stdev.unscaled)
SE_admixture <- select(SE, "ConditionFlu:Admixture", "ConditionNI:Admixture")

mash_data = mash_set_data(as.matrix(betas_admixture), as.matrix(SE_admixture))

############################################
## Step 2: set up the covariance matrices ##
############################################
## use the  “canonical” covariance matries
U.c = cov_canonical(mash_data)
print(names(U.c))

## estimate null correlations and add to data (decreases number of false positive from my test)
V = estimate_null_correlation_simple(mash_data)
data_V = mash_update_data(mash_data, V = V)

## set up the data-driven covariance matrix
m.1by1 = mash_1by1(data_V)
strong = get_significant_results(m.1by1, 0.05)
betas <- as.matrix(betas_admixture)
strong_betas <- betas[strong,]

# obtain initial data_driven covariance matrixes
U.pca = cov_pca(data_V, 2, subset=strong)
print(names(U.pca))
## extreme deconvolution
U.ed = cov_ed(data_V, U.pca, subset=strong)

############################################
########## Step 3: fit the model ###########
############################################
## This step must be peformed using all the tests (or a large random subset),
## because this is where mash learns that many tests are null and corrects for it.

## run model with both canonical and data-driven cov matrices
m.c = mash(data_V, c(U.c, U.ed))
print(get_loglik(m.c),digits = 10)

############################################
### Step 4: Extract Posterior Summaries ####
############################################
# Get effects that are “significant”, which  means they have lfsr less than .10 in at least one condition
thresh <- .10
m.pairwise_PM <- get_pairwise_sharing(m.c, lfsr_thresh=thresh, factor = 0.5)
m.sig <- get_significant_results(m.c, thresh=thresh, sig_fn=get_lfdr)

############################################################
#### Identify condition-specific and shared PopDE genes ####
############################################################

#write this as an output and place in summary directory.
sink(paste0(out_dir, data, "/mash_summary_stats/summary_info.txt"))

print("how many are significant in at least one condition")
print(length(m.sig))
print("how many are significant in just Flu")
length(get_significant_results(m.c, thresh=thresh, sig_fn=get_lfsr, conditions=1))
print("how many are significant in just NI")
length(get_significant_results(m.c, thresh=thresh, sig_fn=get_lfsr, conditions=2))

#finish printing to file.
sink()

## write posterior outs
write.table(get_lfsr(m.c), paste0(out_dir, data, "/mash_results/lfsr_output.txt"), quote = FALSE)
write.table(get_lfdr(m.c), paste0(out_dir, data, "/mash_results/lfdr_output.txt"), quote = FALSE)
write.table(get_pm(m.c), paste0(out_dir, data, "/mash_results/posteriorMeans.txt"), quote = FALSE)
write.table(get_psd(m.c), paste0(out_dir, data, "/mash_results/posteriorStandardDevs.txt"), quote = FALSE)

## save R object
save(m.c, file=paste0(out_dir, data, "/mash_results/mash_output.Rdata"))

## Write results tables with lfdr - local FDR
lfdr <- as.data.frame(get_lfdr(m.c))
lfdr_FLU <- select(lfdr, "ConditionFlu:Admixture")
colnames(lfdr_FLU) <- c("lfdr")

results_FLU_lfdr <- merge(results_FLU, lfdr_FLU, by=0)
rownames(results_FLU_lfdr) <- results_FLU_lfdr$Row.names

## also merge lfsr - local false sign rate (more conservative).
lfsr <- as.data.frame(get_lfsr(m.c))
lfsr_FLU <- select(lfsr, "ConditionFlu:Admixture")
colnames(lfsr_FLU) <- c("lfsr")

results_FLU_lfdr_lfsr <- merge(results_FLU_lfdr, lfsr_FLU, by=0)
rownames(results_FLU_lfdr_lfsr) <- results_FLU_lfdr_lfsr$Row.names


## now for NI
lfdr_NI<- select(lfdr, "ConditionNI:Admixture")
colnames(lfdr_NI) <- c("lfdr")

results_NI_lfdr <- merge(results_NI, lfdr_NI, by=0)
rownames(results_NI_lfdr) <- results_NI_lfdr$Row.names

## also	merge lfsr - local false sign rate (more conservative).
lfsr_NI <- select(lfsr, "ConditionNI:Admixture")
colnames(lfsr_NI) <- c("lfsr")

results_NI_lfdr_lfsr <- merge(results_NI_lfdr, lfsr_NI, by=0)
rownames(results_NI_lfdr_lfsr) <- results_NI_lfdr_lfsr$Row.names

# Write mash results
write.table(results_FLU_lfdr_lfsr, paste0(out_dir, data, "/mash_results/popDE_FLU.txt"))
write.table(results_NI_lfdr_lfsr, paste0(out_dir, data, "/mash_results/popDE_NI.txt"))