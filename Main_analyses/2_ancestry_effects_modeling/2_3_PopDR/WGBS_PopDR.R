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
library(preprocessCore)
library(limma)
library(edgeR)
library(statmod)
library(sva)

data <- "WGBS"
## Use previously produced permutations (TRUE or FALSE)
reproduce <- FALSE

## Set directory structure
folder = "2_ancestry_effects_modeling"
## Set directory 3 steps above script.
setwd('../../../')
system(paste0("mkdir -p Outputs/",folder,"/PopDR_analysis/"))
out_dir <- paste0("Outputs/",folder,"/PopDR_analysis/")

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

#Load DE infection results
DE_infection_results <- read.table(paste0("Outputs/1_DE_infection_modeling/", data, "_results/resultsALL.txt"), header=TRUE,  row.names = 1, sep=",")

#subset to only include those samples for which there is a 1 in the "PopDEset" column (gives you an opportunity to exclude any samples)
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


###########################################################
### Format reads matrices & metadata for PopDR analysis ###
###########################################################

## Get counts
cts <- getMeth(BS.fit, type = "raw")

### Get the condition-specific reads matrices
##subset the NI samples only)
reads_NI=cts[,which(meta_data$Condition=="NI")]
##subset the Flu samples only)
reads_FLU=cts[,which(meta_data$Condition=="Flu")]

### Get condition-specific metadata tables
cols_NI=meta_data[which(meta_data$Condition=="NI"),]
cols_FLU=meta_data[which(meta_data$Condition=="Flu"),]

#this should be 0. It is a check to ensure that all the metadata and read counts matrix have all the same samples.
length(which(rownames(cols_NI)!=colnames(reads_NI)))
length(which(rownames(cols_FLU)!=colnames(reads_FLU)))


##################################################################################
### Build Fold-change matrices & weights & metadata tables from paired samples ###
##################################################################################

#Fold-change function for WGBS cts.
##Notes: FC will be NA if one of the conditions has an NA.
get_FC=function(exp_stim,exp_ref,cols_stim,cols_ref,cond_stim,cond_ref){

  FC=exp_stim

  cont=0

  cols_all=rbind(cols_ref,cols_stim)
  exp_all=cbind(exp_ref,exp_stim)

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
      colnames(FC)[cont]=levels(cols_all$Genotyping_ID)[i]
    }
  }

  FC=FC[,1:cont]
  return(FC)
}

## Get FC matrix from WGBS counts
FC_FLU=get_FC(exp_stim=reads_FLU,exp_ref=reads_NI,cols_stim=cols_FLU,cols_ref=cols_NI,cond_stim="Flu",cond_ref="NI")

## Get metadata for  FC matrix
cols_FLU_FC=meta_data[which(meta_data$Condition=="Flu"),]
cols_FLU_FC=cols_FLU_FC[which(cols_FLU_FC$Genotyping_ID %in% colnames(FC_FLU)),]
rownames(cols_FLU_FC)=cols_FLU_FC$Genotyping_ID

## # should be 0 if everything is concordant
length(which(cols_FLU_FC$Genotyping_ID!=colnames(FC_FLU)))


#################################################################
###  Model PopDR effects of admixture on responses to ligands ###
#################################################################

##  Subset to include only those features which had a q<.20 for Infection DE change
sig_infection_effects <- DE_infection_results[DE_infection_results$qvals <.20, ]
## how many significant features are there?
dim(sig_infection_effects)[1]
## set rownames
sig_infection_effects <- rownames(sig_infection_effects)
## subset FC matrix to only include those that were significant
infection_sig_FC_FLU <- subset(FC_FLU, rownames(FC_FLU) %in% sig_infection_effects)
# check (should be 0)
length(which(rownames(infection_sig_FC_FLU)!=sig_infection_effects))

## Need to also include Age and Batch in model here
DR=function(cols,exp)
{
  design=model.matrix(~Admixture+Age+Batch,data=cols)
  fit <-lmFit(exp,design)
  fit <- eBayes(fit)
  return(fit)
}

DR_FLU=DR(cols_FLU_FC,infection_sig_FC_FLU)

#create a data.frame of the results (and remove automatic p-value adjustment calculation)
results_FLU=topTable(DR_FLU,coef="Admixture",number=nrow(DR_FLU$coefficients))[,1:4]


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

    FLU_rand=DR(cols=cols_FLU_FC_random,exp=infection_sig_FC_FLU)

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