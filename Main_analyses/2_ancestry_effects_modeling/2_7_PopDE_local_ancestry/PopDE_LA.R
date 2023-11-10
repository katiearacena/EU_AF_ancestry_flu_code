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

# Perform 5 permutations
iterations=5

#############################################
## CREATE DIRECTORY STRUCTURE & LOAD DATA  ##
#############################################

## Set directory structure
folder = "2_ancestry_effects_modeling"
## Set directory 3 steps above script.
setwd('../../../')
system(paste0("mkdir -p Outputs/",folder,"/PopDE_local_ancestry/", data, "/"))
out_dir <- paste0("Outputs/",folder,"/PopDE_local_ancestry/", data, "/")

#Set inputs

#Load metadata
meta_data = read.table(paste0("Inputs/metadata/", data, "_metadata.txt"))
#Load read counts (row 1 as rownames)
reads = read.table(paste0("Outputs/", folder,"/batch_corrected_cts_matrices/",data, "_batch.age.corrected.txt"), header=TRUE, row.names=1)
#Load local ancestry estimates
LA_matrix = read.table(paste0("Inputs/LA/", data, "_matrix_for_LA.txt"))
LA_matrix = select(LA_matrix, -c(LA_region))

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

## edit LA calls to have columns for both conditions
LA_geno_NI <- LA_matrix
LA_geno_Flu <- LA_matrix
colnames(LA_geno_NI) <- paste0(colnames(LA_geno_NI), "_NI")
colnames(LA_geno_Flu) <- paste0(colnames(LA_geno_Flu), "_Flu")

LA_matrix_both_conds <- merge(LA_geno_NI, LA_geno_Flu, by.x="feature_NI", by.y="feature_Flu")

## set feature or genes as rownames and remove feature column
if(data=="RNAseq"){
    ## remove those genes that are not in reads (due to using all ensemble id positions)
    LA_matrix_both_conds <- LA_matrix_both_conds[LA_matrix_both_conds$Gene_ID_NI %in% rownames(reads), ]
    rownames(LA_matrix_both_conds) <- LA_matrix_both_conds$Gene_ID_NI
    LA_matrix_both_conds <- select(LA_matrix_both_conds, -c(Gene_ID_NI, Gene_ID_Flu, feature_NI))
    LA_matrix_both_conds <- LA_matrix_both_conds[complete.cases(LA_matrix_both_conds), ]
} else {
    rownames(LA_matrix_both_conds) <- LA_matrix_both_conds$feature_NI
    LA_matrix_both_conds <- select(LA_matrix_both_conds, -c(feature_NI))
}



#####################################################
#### Model DE within each condition (real data) #####
#####################################################

## set function. Compared to the original popDE model "Admixture" is replaced with "LA_call"
DE_WC=function(reads,cols){
  design = model.matrix(~0+ Condition + Condition:LA_call, data = cols)
  # Remove mock columns
  design <- design[, colSums(design != 0) > 0]
  fit <-lmFit(reads,design)
  fit <- eBayes(fit)
  return(list(fit,reads))
}

## for each feature
features <- rownames(reads)

for(j in 1:length(features)){
  if(j%%1000 == 0) print(j)
  feature <- features[j]

  LA_call <- as.data.frame(t(LA_matrix_both_conds[rownames(LA_matrix_both_conds) %in% feature,]))

  if(dim(LA_call)[2]>0){

    colnames(LA_call)[1] <- "LA_call"
    LA_call$indiv_ID <- rownames(LA_call)
    ## should be 0.
    length(which(rownames(LA_call) != LA_call$indiv_ID))
    meta_data_loop <- merge(meta_data, LA_call, by=0)
    meta_data_loop$LA_call <- as.numeric(meta_data_loop$LA_call)
    rownames(meta_data_loop) <- meta_data_loop$Row.names;  meta_data_loop$Row.names <- NULL

    DE_BOTH=DE_WC(reads=reads,cols=meta_data_loop)
    fit_BOTH=DE_BOTH[[1]]
    v_BOTH=DE_BOTH[[2]]

    ## for each feature, store each output separately
    results_NI=topTable(fit_BOTH[rownames(fit_BOTH) %in% feature , ],coef="ConditionNI:LA_call",number=nrow(fit_BOTH$coefficients))[,1:4]
    results_FLU=topTable(fit_BOTH[rownames(fit_BOTH) %in% feature , ],coef="ConditionFlu:LA_call",number=nrow(fit_BOTH$coefficients))[,1:4]
        
  } else if (dim(LA_call)[2]==0) {

    results_NI <- data.frame(NA, NA, NA, NA)
    results_FLU <- data.frame(NA, NA, NA, NA)
    colnames(results_NI) <- c("logFC", "AveExpr", "t", "P.Value")
    colnames(results_FLU) <- c("logFC", "AveExpr", "t", "P.Value")
    rownames(results_NI) <- feature
    rownames(results_FLU) <- feature
  }

  if(j == 1){
    NI_outs <- results_NI
    FLU_outs <- results_FLU
  }else{
    NI_outs <- rbind(NI_outs, results_NI)
    FLU_outs <- rbind(FLU_outs, results_FLU)
  }
}

#create the directory structure for the output files.
system(paste0("mkdir -p ", out_dir, "/results"))
system(paste0("mkdir -p ", out_dir, "/summary_stats/"))

#Write outputs
write.table(NI_outs, paste0(out_dir, "/results/popDE_NI.LA.real.txt"))
write.table(FLU_outs, paste0(out_dir, "/results/popDE_FLU.LA.real.txt"))


##############################################################
#### Do permutation tests to define null model for PopDE #####
##############################################################

if(reproduce)
{
  ### to reproduce results previously generated load the permuted data.
  shuffled_pvals_adm_NI=read.table(paste0(out_dir,"/permuted_pvalues/popDE_NI.txt"))
  shuffled_pvals_adm_FLU=read.table(paste0(out_dir,"/permuted_pvalues/popDE_FLU.txt"))
  ##if reproduce is set to TRUE this part below will be skipped. If reproduce = FALSE then you the code below will permute the data.

}else{

  for(iter in 1:iterations)
  {
    print(iter)

    for(j in 1:length(features)){
      if(j%%1000 == 0) print(paste0("perm ",iter,": ", j))
      feature <- features[j]
      
      LA_call <- as.data.frame(t(LA_matrix_both_conds[rownames(LA_matrix_both_conds) %in% feature,]))

      ## shuffle labels to permute the data
      LA_genotype_to_merge_SHUFFLED= LA_call
      LA_genotype_to_merge_SHUFFLED[,1] = sample(LA_genotype_to_merge_SHUFFLED[,1])


      if(dim(LA_genotype_to_merge_SHUFFLED)[2]>0){
        colnames(LA_genotype_to_merge_SHUFFLED)[1] <- "LA_call"

        ## add these genotype into the metadata
        meta_data_loop <- meta_data
        meta_data_loop <- merge(meta_data_loop, LA_genotype_to_merge_SHUFFLED, by=0)
        meta_data_loop$LA_call <- as.numeric(meta_data_loop$LA_call)
        rownames(meta_data_loop) <- meta_data_loop$Row.names;  meta_data_loop$Row.names <- NULL

        DE_BOTH=DE_WC(reads=reads,cols=meta_data_loop)
        fit_BOTH=DE_BOTH[[1]]
        v_BOTH=DE_BOTH[[2]]

        ## for each feature, store each output separately
        perm_NI=topTable(fit_BOTH[rownames(fit_BOTH) %in% feature , ],coef="ConditionNI:LA_call",number=nrow(fit_BOTH$coefficients))[,1:4]
        perm_FLU=topTable(fit_BOTH[rownames(fit_BOTH) %in% feature , ],coef="ConditionFlu:LA_call",number=nrow(fit_BOTH$coefficients))[,1:4]
      
      } else if (dim(LA_genotype_to_merge_SHUFFLED)[2]==0) {
        perm_NI <- data.frame(NA, NA, NA, NA)
        perm_FLU <- data.frame(NA, NA, NA, NA)
        colnames(perm_NI) <- c("logFC", "AveExpr", "t", "P.Value")
        colnames(perm_FLU) <- c("logFC", "AveExpr", "t", "P.Value")
        rownames(perm_NI) <- feature
        rownames(perm_FLU) <- feature
      }

      if(j == 1){
        NI_perm_outs <- perm_NI
        FLU_perm_outs <- perm_FLU
      }else{
        NI_perm_outs <- rbind(NI_perm_outs, perm_NI)
        FLU_perm_outs <- rbind(FLU_perm_outs, perm_FLU)
      }
    }

    permuted_adm_FLU=FLU_perm_outs
    permuted_adm_NI=NI_perm_outs

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
  system(paste0("mkdir -p ", out_dir, "/permuted_pvalues"))
  #Writing permutation tests p-values tables
  write.table(shuffled_pvals_adm_FLU,paste0(out_dir, "/permuted_pvalues/popDE.LA_FLU.txt"))
  write.table(shuffled_pvals_adm_NI,paste0(out_dir, "/permuted_pvalues/popDE.LA_NI.txt"))

}

##################################################
#### Write results and perform FDR correction ####
##################################################

## Remove NAs (will be same between real data and permutations)
FLU_outs <- FLU_outs[complete.cases(FLU_outs), ]
NI_outs <- NI_outs[complete.cases(NI_outs), ]

shuffled_pvals_adm_FLU <- shuffled_pvals_adm_FLU[complete.cases(shuffled_pvals_adm_FLU), ]
shuffled_pvals_adm_NI <- shuffled_pvals_adm_NI[complete.cases(shuffled_pvals_adm_NI), ]

dim(FLU_outs)[1] == dim(shuffled_pvals_adm_FLU)[1]
dim(NI_outs)[1] == dim(shuffled_pvals_adm_NI)[1]

##calculate BH adjusted p-values.
BH_Flu_pvals <- p.adjust(FLU_outs$P.Value, method= "BH", n=length(FLU_outs$P.Value))
BH_NI_pvals <- p.adjust(NI_outs$P.Value, method= "BH", n=length(NI_outs$P.Value))

##Bind BH adjusted p-values to results table.
results_FLU_BH=cbind(FLU_outs, BH_Flu_pvals)
results_NI_BH=cbind(NI_outs, BH_NI_pvals)

## Calculate qvalue
q_Flu_pvals <- empPvals(stat=-log10(FLU_outs$P.Value), stat0=-log10(as.matrix(shuffled_pvals_adm_FLU)), pool = T)
qvals_Flu <- qvalue(p=q_Flu_pvals)$qvalue
q_NI_pvals <- empPvals(stat=-log10(NI_outs$P.Value), stat0=-log10(as.matrix(shuffled_pvals_adm_NI)), pool = T)
qvals_NI <- qvalue(p=q_NI_pvals)$qvalue

## Bind qvalue to same tables.
results_FLU_BH_q=cbind(results_FLU_BH, qvals_Flu)
results_NI_BH_q=cbind(results_NI_BH, qvals_NI)

## Write results tables
write.table(results_FLU_BH_q, paste0(out_dir, "/results/popDE.LA_FLU.txt"))
write.table(results_NI_BH_q, paste0(out_dir, "/results/popDE.LA_NI.txt"))

###################################
#### Identify sig. PopDE genes ####
###################################

#write this as an output and place in summary directory.
sink(paste0(out_dir, "/summary_stats/summary_info.txt"))

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
system(paste0("mkdir -p ", out_dir, "/mash_results"))
system(paste0("mkdir -p ", out_dir, "/mash_summary_stats"))


############################################
######### Step 1: Read in the data #########
############################################
## To run mash you need data consisting of a matrix of effects (Bhat) and a matrix of standard errors (Shat),
## for 𝐽 effects (rows) in 𝑅 conditions (columns).

##extract betas for mash
Flu_res_betas <- select(results_FLU_BH_q, logFC)
NI_res_betas <- select(results_NI_BH_q, logFC)
betas_LA <- merge(Flu_res_betas, NI_res_betas, by=0)
rownames(betas_LA) <- betas_LA$Row.names; betas_LA$Row.names <- NULL
colnames(betas_LA) <- c("betas_Flu", "betas_NI")

## Get standard errors of betas 
results_NI_BH_q$beta_se = results_NI_BH_q$logFC/results_NI_BH_q$t
results_FLU_BH_q$beta_se = results_FLU_BH_q$logFC/results_FLU_BH_q$t

Flu_res_ses <- select(results_FLU_BH_q, beta_se)
NI_res_ses <- select(results_NI_BH_q, beta_se)
ses_LA <- merge(Flu_res_ses, NI_res_ses, by=0)
rownames(ses_LA) <- ses_LA$Row.names; ses_LA$Row.names <- NULL
colnames(ses_LA) <- c("ses_Flu", "ses_NI")

mash_data = mash_set_data(as.matrix(betas_LA), as.matrix(ses_LA))

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
betas <- as.matrix(betas_LA)
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
sink(paste0(out_dir, "/mash_summary_stats/summary_info.txt"))

print("how many are significant in at least one condition")
print(length(m.sig))
print("how many are significant in just Flu")
length(get_significant_results(m.c, thresh=thresh, sig_fn=get_lfsr, conditions=1))
print("how many are significant in just NI")
length(get_significant_results(m.c, thresh=thresh, sig_fn=get_lfsr, conditions=2))

#finish printing to file.
sink()

## write posterior outs
write.table(get_lfsr(m.c), paste0(out_dir, "/mash_results/lfsr_output.txt"), quote = FALSE)
write.table(get_lfdr(m.c), paste0(out_dir, "/mash_results/lfdr_output.txt"), quote = FALSE)
write.table(get_pm(m.c), paste0(out_dir, "/mash_results/posteriorMeans.txt"), quote = FALSE)
write.table(get_psd(m.c), paste0(out_dir, "/mash_results/posteriorStandardDevs.txt"), quote = FALSE)

## save R object
save(m.c, file=paste0(out_dir, "/mash_results/mash_output.Rdata"))

## Write results tables with lfdr - local FDR
lfdr <- as.data.frame(get_lfdr(m.c))
lfdr_FLU <- select(lfdr, "betas_Flu")
colnames(lfdr_FLU) <- c("lfdr")

results_FLU_lfdr <- merge(results_FLU, lfdr_FLU, by=0)
rownames(results_FLU_lfdr) <- results_FLU_lfdr$Row.names

## also merge lfsr - local false sign rate (more conservative).
lfsr <- as.data.frame(get_lfsr(m.c))
lfsr_FLU <- select(lfsr, "betas_Flu")
colnames(lfsr_FLU) <- c("lfsr")

results_FLU_lfdr_lfsr <- merge(results_FLU_lfdr, lfsr_FLU, by=0)
rownames(results_FLU_lfdr_lfsr) <- results_FLU_lfdr_lfsr$Row.names


## now for NI
lfdr_NI<- select(lfdr, "betas_NI")
colnames(lfdr_NI) <- c("lfdr")

results_NI_lfdr <- merge(results_NI, lfdr_NI, by=0)
rownames(results_NI_lfdr) <- results_NI_lfdr$Row.names

## also	merge lfsr - local false sign rate (more conservative).
lfsr_NI <- select(lfsr, "betas_NI")
colnames(lfsr_NI) <- c("lfsr")

results_NI_lfdr_lfsr <- merge(results_NI_lfdr, lfsr_NI, by=0)
rownames(results_NI_lfdr_lfsr) <- results_NI_lfdr_lfsr$Row.names

# Write mash results
write.table(results_FLU_lfdr_lfsr, paste0(out_dir, "/mash_results/popDE.LA_FLU.txt"))
write.table(results_NI_lfdr_lfsr, paste0(out_dir, "/mash_results/popDE.LA_NI.txt"))
