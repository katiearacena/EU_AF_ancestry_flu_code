#load libraries
library(preprocessCore)
library(ggplot2)
library(limma)
library(edgeR)
library(statmod)
library(sva)
library(dplyr)
library(qvalue)


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
system(paste0("mkdir -p Outputs/",folder,"/PopDR_analysis/"))
out_dir <- paste0("Outputs/",folder,"/PopDR_analysis/")

#Set inputs

#Load metadata
meta_data = read.table(paste0("Inputs/metadata/", data, "_metadata.txt"))
#Load read counts (row 1 as rownames)
reads_NI = read.table(paste0("Outputs/", folder,"/batch_corrected_cts_matrices/",data, "_batch.age.corrected_NI.samples.txt"), header=TRUE, row.names=1)
reads_FLU = read.table(paste0("Outputs/", folder,"/batch_corrected_cts_matrices/",data, "_batch.age.corrected_FLU.samples.txt"), header=TRUE, row.names=1)

#Load DE infection results
DE_infection_results <- read.table(paste0("Outputs/1_DE_infection_modeling/", data, "_results/resultsALL.txt"), header=TRUE,  row.names = 1, sep=",")


########################################
### Format reads matrices & metadata ###
########################################

#subset to only include those samples for which there is a 1 in the "PopDEset" column
meta_data=meta_data[which(meta_data$PopDE_set==1),]
#ensure that the Sample_ID are the rownames of the cols dataframe
rownames(meta_data)=meta_data$Sample
#ensure that the cols dataframe is ordered alphabetically (based on sample_id)
meta_data=meta_data[order(rownames(meta_data)),]

#remove data type string from chipseq metadata rownames and Sample ID
rownames(meta_data)<- sub("_H3K27ac", "", rownames(meta_data))
rownames(meta_data)<- sub("_H3K4me1", "", rownames(meta_data))
rownames(meta_data)<- sub("_H3K27me3", "", rownames(meta_data))
rownames(meta_data)<- sub("_H3K4me3", "", rownames(meta_data))

meta_data$Sample<- sub("_H3K27ac", "", meta_data$Sample)
meta_data$Sample<- sub("_H3K4me1", "", meta_data$Sample)
meta_data$Sample<- sub("_H3K27me3", "", meta_data$Sample)
meta_data$Sample<- sub("_H3K4me3", "", meta_data$Sample)

#make sure that all the columns of the meta data are also in the read counts matrix.
reads_NI=reads_NI[,which(colnames(reads_NI) %in% rownames(meta_data))]
reads_FLU=reads_FLU[,which(colnames(reads_FLU) %in% rownames(meta_data))]

#order the reads count matrix.
reads_NI=reads_NI[,order(colnames(reads_NI))]
reads_FLU=reads_FLU[,order(colnames(reads_FLU))]

### Get condition-specific metadata tables
cols_NI=meta_data[which(meta_data$Condition=="NI"),]
cols_FLU=meta_data[which(meta_data$Condition=="Flu"),]

#this should be 0. It is a check to ensure that all the metadata and read counts matrix have all the same samples.
length(which(rownames(cols_NI)!=colnames(reads_NI)))
length(which(rownames(cols_FLU)!=colnames(reads_FLU)))

# Load weights and order
weights_NI=read.table(paste0("Outputs/", folder,"/batch_corrected_cts_matrices/",data, "_weights_for_PopDR_NI.samples.txt"), header=TRUE)
weights_FLU=read.table(paste0("Outputs/", folder,"/batch_corrected_cts_matrices/",data, "_weights_for_PopDR_FLU.samples.txt"), header=TRUE)

weights_NI=weights_NI[,order(colnames(reads_NI))]
weights_FLU=weights_FLU[,order(colnames(reads_FLU))]

####################################################################################
#### Build Fold-change matrices & weights & metadata tables from paired samples ####
####################################################################################

# Set fold-change function
get_FC=function(exp_stim,exp_ref,weights_stim,weights_ref,cols_stim,cols_ref,cond_stim,cond_ref){

  FC=exp_stim
  weights_FC=FC

  cont=0

  cols_all=rbind(cols_ref,cols_stim)
  exp_all=cbind(exp_ref,exp_stim)
  weights_all=cbind(weights_ref,weights_stim)

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
      weights_FC[,cont]=1/(1/weights_all[,stim_index]+1/weights_all[,ref_index])

      colnames(FC)[cont]=levels(cols_all$Genotyping_ID)[i]
      colnames(weights_FC)[cont]=levels(cols_all$Genotyping_ID)[i]
    }
  }

  FC=FC[,1:cont]
  weights_FC=weights_FC[,1:cont]

  return(list(FC,weights_FC))
}

## Apply function 
FC_FLU=get_FC(exp_stim=reads_FLU,exp_ref=reads_NI,weights_stim=weights_FLU,weights_ref=weights_NI,cols_stim=cols_FLU,cols_ref=cols_NI,cond_stim="Flu",cond_ref="NI")

## Extract FC matrix and weights
weights_FLU=FC_FLU[[2]]
FC_FLU=FC_FLU[[1]]

cols_FLU_FC=meta_data[which(meta_data$Condition=="Flu"),]
cols_FLU_FC=cols_FLU_FC[which(cols_FLU_FC$Genotyping_ID %in% colnames(FC_FLU)),]
rownames(cols_FLU_FC)=cols_FLU_FC$Genotyping_ID

# should be 0 if everything is concordant.
length(which(cols_FLU_FC$Genotyping_ID!=colnames(FC_FLU)))
length(which(cols_FLU_FC$Genotyping_ID!=colnames(weights_FLU)))


##################################################################
#### Model PopDR effects of admixture on responses to ligands ####
##################################################################

##  Subset to include only those features which had a q<.10 for Infection DE change
sig_infection_effects <- DE_infection_results[DE_infection_results$qvals <.10, ]
## how many significant features are there?
dim(sig_infection_effects)[1]
## set rownames
sig_infection_effects <- rownames(sig_infection_effects)
## subset FC matrix to only include those that were significant
infection_sig_FC_FLU <- subset(FC_FLU, rownames(FC_FLU) %in% sig_infection_effects)
infection_sig_weights_FLU <- subset(weights_FLU, rownames(weights_FLU) %in% sig_infection_effects)
# check (should be 0)
length(which(rownames(infection_sig_FC_FLU)!=sig_infection_effects))
length(which(rownames(infection_sig_weights_FLU)!=sig_infection_effects))

## Set PopDR model
DR=function(cols,exp,weights)
{
  design=model.matrix(~Admixture,data=cols)
  fit <-lmFit(exp,weights=weights,design)
  fit <- eBayes(fit)
  return(fit)
}

## model PopDR effects
DR_FLU=DR(cols_FLU_FC,infection_sig_FC_FLU,infection_sig_weights_FLU)


#################################################################################
#### Randomize admixture in permutation tests to define null model for PopDR ####
#################################################################################

## permute 1000 times
iterations=1000

if(reproduce)
{
  ### to reproduce results previously generated load the permuted data.
  shuffled_pvals_adm_FLU=read.table(paste0(out_dir, data,"/permuted_pvalues/popDR_FLU.txt"))
}else{

  cols_FLU_FC_random=cols_FLU_FC
  for(iter in 1:iterations)
  {
    print(iter)
    ## permute admixture 
    cols_FLU_FC_random$Admixture=sample(cols_FLU_FC_random$Admixture)

    FLU_rand=DR(cols=cols_FLU_FC_random,exp=infection_sig_FC_FLU,weights=infection_sig_weights_FLU)

    permuted_adm_FLU=topTable(FLU_rand, coef="Admixture", sort.by="none", adjust="BH",n=nrow(FLU_rand$coefficients))

    if(iter==1)
    {
      shuffled_pvals_adm_FLU <-data.frame(x=permuted_adm_FLU$P.Value)

      rownames(shuffled_pvals_adm_FLU)=rownames(permuted_adm_FLU)
    } else {
      shuffled_pvals_adm_FLU <- cbind(shuffled_pvals_adm_FLU,x=permuted_adm_FLU$P.Value)

    }
  }

  system(paste0("mkdir -p ", out_dir, data, "/permuted_pvalues"))
  #Writing permutation tests p-values tables
  write.table(shuffled_pvals_adm_FLU,paste0(out_dir, data,"/permuted_pvalues/popDR_FLU.txt"))
}

#create the directory structure for the output files.
system(paste0("mkdir -p ", out_dir, data, "/results"))
system(paste0("mkdir -p ", out_dir, data, "/summary_stats/"))

#create a data.frame of the results (and remove automatic p-value adjustment calculation)
results_FLU=topTable(DR_FLU,coef="Admixture",number=nrow(DR_FLU$coefficients))[,1:4]

##################################################
#### Write results and perform FDR correction ####
##################################################

##calculate BH adjusted p-values.
BH_FLU_pvals <- p.adjust(results_FLU$P.Value, method= "BH", n=length(results_FLU$P.Value))

## Calculate qvalue
q_FLU_pvals <- empPvals(stat=-log10(results_FLU$P.Value), stat0=-log10(as.matrix(shuffled_pvals_adm_FLU)), pool = T)
qvals_FLU <- qvalue(p=q_FLU_pvals)$qvalue

## Bind BH and qvalue to results table.
results_FLU_BH_q=cbind(results_FLU, BH_FLU_pvals, qvals_FLU)

## Write outputs.
write.table(x=results_FLU_BH_q, file= paste0(out_dir,data,"/results/results_FLU.txt"))
## Write original FC matrix for plotting examples.
write.table(x=FC_FLU, file= paste0(out_dir,data,"/results/FC_for_plotting.txt"))

#what are the PopDE genes in each condition using the different FDR thresholds?
#write this as an output and place in summary  directory
sink(paste0(out_dir, data, "/summary_stats/summary_info.txt"))

print("what are the PopDR genes in each condition using the different FDR thresholds?")

print("sig hits with BH FDR <= 0.20")
print(summary(results_FLU_BH_q$BH_FLU_pvals <= 0.20))
print("sig hits with q <= 0.20")
print(summary(results_FLU_BH_q$qvals_FLU <= 0.20))

#finish printing to file.
sink()