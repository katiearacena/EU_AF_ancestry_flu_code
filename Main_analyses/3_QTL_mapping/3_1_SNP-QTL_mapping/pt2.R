#The following script is modified from Harrison et al. 2019.


##############################################################
#### 1. Load dependencies ####################################
##############################################################

## Datatype
condition <- "NI"
## Declare number of expression PCs to regress out.
expPCs_reg <-"1"
## Declare directory for temp files. Include PCs reg and condition.
temp_dir <- "1expPC_NI"
## Declare output directory.
res_dir <- "1expPCs"


## Set number of random permutations (for empiric fdr corrections).
## Normally, this is set to 10 iterations.
iterations=10

## Load required libraries.
library(MatrixEQTL)
library(edgeR)
library(limma)
library(gdsfmt)
library(SNPRelate)
library(qvalue)
library(gplots)
library(viridis)
library(ggplot2)
library(corrplot)
library(rmeta)
library(bsseq)


data <- "WGBS"
# Perform 10 permutations to use for FDR correction.
iterations <- 10

#############################################################
## CREATE DIRECTORY STRUCTURE, SET FUNCTIONS, & LOAD DATA  ##
#############################################################

folder = "3_QTL_mapping"
## Set directory 3 steps above script.
setwd('../../../')

## Create directory structure to save outputs.

out_dir <- paste0("Outputs/",folder,"/SNP-QTL_mapping/", data, "/", condition, "/", res_dir, "/")

permuted_pvalues_folder=paste0(out_dir,"/raw_results/")

#######################################################################################################
## Write Best SNP-gene associations files for Delta-PVE analyses & prepare files for FDR corrections ##
#######################################################################################################

## Select the top cis-SNP for each gene in true and permuted files
for(iteration in 0:iterations)
{   print(iteration)
    if(iteration==0){
        event=read.table(paste0(out_dir,"/raw_results/result_original.txt"),header=TRUE)
        event.sort<-event[order(event[,2],event[,4]),]
        event.bestQTL<-event.sort[!duplicated(event.sort$gene),]
        event.bestQTL<-event.bestQTL[order(event.bestQTL[,4]),]
    }else{
        event=read.table(paste0(permuted_pvalues_folder,"/result_permuted_",iteration,".txt"),header=TRUE)
        event.sort<-event[order(event[,2],event[,4]),]
        event.bestQTL<-event.sort[!duplicated(event.sort$gene),]
        event.bestQTL<-event.bestQTL[order(event.bestQTL[,5]),]
    }
    if(iteration==0){
        original_best_EQTL=event.bestQTL
    }else{
        if(iteration==1)
        {
            permuted1_best_EQTL=event.bestQTL
            #only pull out the p-value.
            Permutation_Input = event.bestQTL[4]
        }else{
            Permutation_Input=cbind(Permutation_Input,event.bestQTL[4])}
    }
}

print("finished loading results and permutations")

#########################################
## Use qvalue package to calculate FDR ##
#########################################

## Calculate qvalues using the best_SNPs files.
emp_pvalues <- empPvals(stat=-log10(original_best_EQTL[,4]), stat0=-log10(as.matrix(Permutation_Input)), pool = T)
qvalues <- qvalue(p=emp_pvalues)$qvalue

print("calculated qvalues")

## Bind qvalue output with results_best_SNPs file
results_best_SNPs_with_qval <- cbind(original_best_EQTL, qvalues)
## Write table with this column added.
write.table(x=results_best_SNPs_with_qval, file= paste0(out_dir,"results_best_SNPs_with_qval.txt"), quote=FALSE)

## Save permuted best SNPs for 1 permutation
write.table(permuted1_best_EQTL, file = paste0(out_dir,"permuted1_best_SNPs.txt"), quote=FALSE)


####################################################
## Create qqplot using p-values of best SNPs only ##
####################################################

## Make and save best SNPs qqplot.
png(filename =  paste0(out_dir,"best_SNPs_qqplot.png"))
qqplot(-log10(permuted1_best_EQTL[,4]), -log10(results_best_SNPs_with_qval[,4]), ylim=  c(0, 40), xlim= c(0, 40), main ="Best SNPs qqplot", xlab = "-log10(theoretical p-values)", ylab = "-log10(observed p-values)")
abline(c(0,1),col="red")
