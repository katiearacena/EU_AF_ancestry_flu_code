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
## perm number
perm_number <- args[2]

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


####################################################
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

##############################################################
#### Do permutation tests to define null model for PopDE #####
#### Save each test and load seperately later            #####
##############################################################

for(j in 1:length(features)){
  if(j%%1000 == 0) print(paste0("perm ",perm_number,": ", j))
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

## save these outputs
permuted_adm_FLU=FLU_perm_outs
permuted_adm_NI=NI_perm_outs

system(paste0("mkdir -p ", out_dir, "/permuted_pvalues"))

#Writing permutation tests p-values tables
write.table(permuted_adm_FLU,paste0(out_dir, "/permuted_pvalues/perm", perm_number, "_popDE.LA_FLU.txt"))
write.table(permuted_adm_NI,paste0(out_dir, "/permuted_pvalues/perm", perm_number, "_popDE.LA_NI.txt"))