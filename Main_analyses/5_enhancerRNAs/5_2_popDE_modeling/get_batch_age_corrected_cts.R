#########################################################################################################
## Voom normalize, remove batch with ComBat, and regress out age effects for count datasets (not WGBS) ##
## This corrected expression matrix is be used for PopDE, PopDR and MatrixeQTL analysis                ##
#########################################################################################################

## Load libraries
library(preprocessCore)
library(limma)
library(edgeR)
library(statmod)
library(sva)
library(reshape2)
library(dplyr)
library(stringr)

## Set directory structure
folder = "5_enhancerRNAs"
## Set directory 3 steps above script.
setwd('../../../')
system(paste0("mkdir -p Outputs/",folder,"/batch_corrected_cts_matrices/"))
outs_dir <- paste0("Outputs/",folder,"/batch_corrected_cts_matrices/")

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


#########################
### Load input files  ###
#########################

#Load metadata
meta_data = read.table(paste0("Inputs/metadata/RNAseq_metadata.txt"))
#Load FILTERED eRNA counts (row 1 as rownames)
reads = read.table(paste0("Outputs/5_enhancerRNAs/DE_infection_effects/eRNA_filtered_reads.txt"), header=TRUE, row.names=1)
## load RNAseq reads
RNAseq_reads = read.table(paste0("Inputs/counts_matrices/RNAseq_filtered.counts.txt"), header=TRUE, row.names=1)


########################################
### Format reads matrices & metadata ###
########################################

#subset to only include those samples for which there is a 1 in the "PopDEset" & "EQTL_set" column.
meta_data=meta_data[which(meta_data$PopDE_set==1 | meta_data$EQTL_set==1),]
#ensure that the Sample_ID are the rownames of the meta_data dataframe
rownames(meta_data)=meta_data$Sample
#ensure that the meta_data dataframe is ordered alphabetically (based on sample_id)
meta_data=meta_data[order(rownames(meta_data)),]
#remove data type string from chipseq metadata

#make sure that all the columns of the meta data are also in the read counts matrix.
reads=reads[,which(colnames(reads) %in% rownames(meta_data))]
RNAseq_reads=RNAseq_reads[,which(colnames(RNAseq_reads) %in% rownames(meta_data))]

#order the reads count matrix.
reads=reads[,order(colnames(reads))]
RNAseq_reads=RNAseq_reads[,order(colnames(RNAseq_reads))]

#this should be 0. It is a check to ensure that all the metadata and read counts matrix have all the same samples.
length(which(rownames(meta_data)!=colnames(reads)))
length(which(rownames(meta_data)!=colnames(RNAseq_reads)))

# Are metadata and reads the same dimensions?
dim(meta_data)[1]==dim(reads)[2]
dim(meta_data)[1]==dim(RNAseq_reads)[2]

# mean center age
meta_data=pretty_up_cols(meta_data)

## bind eRNA reads and RNAseq reads
## make sure colnames are in same order and match
length(which(colnames(RNAseq_reads)!=colnames(reads)))
gene_eRNA_reads <- rbind(reads, RNAseq_reads)
dim(gene_eRNA_reads)

### Get condition-specific reads matrices
##subset the NI samples only)
reads_NI=gene_eRNA_reads[,which(meta_data$Condition=="NI")]
##subset the Flu samples only)
reads_FLU=gene_eRNA_reads[,which(meta_data$Condition=="Flu")]

### Get condition-specific metadata tables
meta_data_NI=refactorize(meta_data[which(meta_data$Condition=="NI"),])
meta_data_FLU=refactorize(meta_data[which(meta_data$Condition=="Flu"),])


######################################################
### Voom normalize reads, remove batch with Combat ###
######################################################

COR=function(reads,meta_data){

  dge <- DGEList(counts=reads)
  dge <- calcNormFactors(dge)
  design=model.matrix(~Age+Admixture,data=meta_data)
  v <- voom(dge,design,plot=FALSE)
  v_combat = ComBat(dat=as.matrix(v$E), batch=meta_data$Batch, mod=design, par.prior=TRUE)
  v$E=v_combat
  fit <-lmFit(v,design)
  fit <- eBayes(fit)
  return(list(fit,v))
}

#Correct for batch within each condition.
COR_NI=COR(reads=reads_NI,meta_data=meta_data_NI)
fit_NI=COR_NI[[1]]
v_NI=COR_NI[[2]]

COR_FLU=COR(reads=reads_FLU,meta_data=meta_data_FLU)
fit_FLU=COR_FLU[[1]]
v_FLU=COR_FLU[[2]]

## Regress out effects of non-admixture covariates within condition

correct_exp=function(v,fit){
  corrected_exp=v$E
  ## regress out second column in the design matrix which is age
  ## the first is intercept and third is admixture
  indexes=c(2)

  for(i in indexes){
    corrected_exp <- corrected_exp - fit$coefficients[,i]%*%t(fit$design[,i])
    return(corrected_exp)
  }
}

corrected_exp_NI=correct_exp(v_NI,fit_NI)
corrected_exp_FLU=correct_exp(v_FLU,fit_FLU)

## bind corrected read counts back together
all_corrected_exp <- cbind(corrected_exp_NI, corrected_exp_FLU)

## subset on only eRNA reads
all_corrected_exp <- all_corrected_exp[grep("chr", rownames(all_corrected_exp)), ]
corrected_exp_FLU <- corrected_exp_FLU[grep("chr", rownames(corrected_exp_FLU)), ]
corrected_exp_NI <- corrected_exp_NI[grep("chr", rownames(corrected_exp_NI)), ]
v_FLU <- v_FLU[grep("chr", rownames(v_FLU)), ]
v_NI <- v_NI[grep("chr", rownames(v_NI)), ]


##############################################
#### Save outputs in appropriate directory ###
##############################################

# Write all corrected expression
write.table(x = all_corrected_exp, file=paste0(outs_dir,"eRNA_batch.age.corrected.txt"), sep="\t")
# Write flu corrected expression
write.table(x = corrected_exp_FLU, file=paste0(outs_dir, "eRNA_batch.age.corrected_FLU.samples.txt"), sep="\t")
# Write NI corrected expression
write.table(x = corrected_exp_NI, file=paste0(outs_dir, "eRNA_batch.age.corrected_NI.samples.txt"), sep="\t")

# Write weights
colnames(v_FLU) <- colnames(corrected_exp_FLU)
write.table(x = v_FLU$weights, file=paste0(outs_dir, "eRNA_weights_for_PopDR_FLU.samples.txt"), sep="\t")

colnames(v_NI) <- colnames(corrected_exp_NI)
write.table(x = v_NI$weights, file=paste0(outs_dir, "eRNA_weights_for_PopDR_NI.samples.txt"), sep="\t")

