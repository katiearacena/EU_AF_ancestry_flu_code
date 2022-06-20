# Load libraries
library(limma)
library(edgeR)
library(ggplot2)
library(cowplot)
library(reticulate)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(gridExtra)
library(ggfortify)
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

#############################################
## CREATE DIRECTORY STRUCTURE & LOAD DATA  ##
#############################################

## Set directory structure
folder = "1_DE_infection_modeling"
## Set directory 3 steps above script.
setwd('../../../')
system(paste0("mkdir -p Outputs/",folder,"/", data, "_residuals"))
system(paste0("mkdir -p Outputs/", folder,"/", data,"_results"))
system(paste0("mkdir -p Outputs/", folder,"/", data,"_voom"))

## Load metadata file
meta_data = read.table(paste0("Inputs/metadata/", data, "_metadata.txt"))
## Remove Mock and IPSC samples which are not used in analysis.
meta_data <- meta_data[which(!meta_data$Condition=="Mock"),]
meta_data <- meta_data[!grepl("EU122", meta_data$Genotyping_ID),]
#remove data type string from chipseq metadata
rownames(meta_data)<- sub("_H3K27ac", "", rownames(meta_data))
rownames(meta_data)<- sub("_H3K4me1", "", rownames(meta_data))
rownames(meta_data)<- sub("_H3K27me3", "", rownames(meta_data))
rownames(meta_data)<- sub("_H3K4me3", "", rownames(meta_data))

#Load read counts (row 1 as rownames)
reads = read.table(paste0("Inputs/counts_matrices/", data, "_filtered.counts.txt"), header=TRUE, row.names=1)

## Set output directories
voom_dir <- paste0("Outputs/", folder,"/", data,"_voom")
results_dir <- paste0("Outputs/", folder,"/", data,"_results")
residuals_dir <- paste0( "Outputs/", folder, "/", data, "_residuals")

#################################
## PREPARE FILES FOR ANALYSIS  ##
#################################

## Factorize certain columns from meta data
meta_data$Batch <- as.factor(meta_data$Batch)
meta_data$Genotyping_ID <- as.factor(meta_data$Genotyping_ID)
## Set NI as first level so that it is used as the reference.
meta_data$Condition = factor(meta_data$Condition, levels=c("NI","Flu"))

## Set original names in case reordering is necessary.
reorder_names <- rownames(meta_data)
length(reorder_names)

## Subset reads to include only the samples in metadata file.
if(length(reorder_names) == dim(reads)[2]){
  reads <- reads[reorder_names]
## If the values do not match, remove those rows from the reads counts file (this will remove Mock samples)
}else{
  correct_names <- colnames(reads)
  meta_data <- meta_data[rownames(meta_data) %in% correct_names,]
  reorder_names <- rownames(meta_data)
  reads <- reads[reorder_names]
}

## Check to make sure reads and meta match 
length(reorder_names) == dim(reads)[2]


#################################
## VOOM NORMALIZE READ COUNTS  ##
#################################

## DGEList simple takes reads and converts to a DGE object format which is needed in downstream functions.
dge <- DGEList(counts = reads)
## CalcNormFactors does normalization by library size (default is TMM)
dge <- calcNormFactors(dge)
## Set design matrix
design = model.matrix(~0 + Batch, data = meta_data)
## Remove any columns that are all 0s
design <- design[, colSums(design != 0) > 0]

## Plot voomed data
pdf(paste0(voom_dir,"/voom_plot.pdf"))
v <- voom(dge, design, plot = TRUE)
dev.off()

## Apply mean variance weights and design matrix to phenotype data.
fit <- lmFit(v, design)
fit <- eBayes(fit)


#################################
##### GET & SAVE RESIDUALS ######
#################################

## Get residuals to regress out batch effect
residuals <- residuals.MArrayLM(object = fit, v)
## residuals should be the same dim as the expression matrix
length(v$E) == length(residuals)

## Calculate average batch effect for each gene
avg_batch_effect <- rowMeans(fit$coefficients)

## Add the average batch effect back into residuals
corrected_expression <- apply(residuals,2,function(x){x + avg_batch_effect})
weights <- v$weights
colnames(weights) <- colnames(corrected_expression)
rownames(weights) <- rownames(corrected_expression)

## Write batch-corrected expression and weights
write.table(corrected_expression, paste0(residuals_dir,"/corrected_expression.txt"), quote = FALSE, sep = ",")
write.table(weights, paste0(residuals_dir,"/weights.txt"), quote = FALSE, sep = ",")


#################################
##### MODEL INFECTION DE  #######
#################################

## Check that corrected_expression and meta_data are concordant
length(which(colnames(corrected_expression)!=rownames(meta_data)))

## Model infection differential expression
design = model.matrix(~0 + Genotyping_ID + Condition, data = meta_data)
design <- design[, colSums(design != 0) > 0]

## Number of pairs used in DE analysis (-1 for condition column)
dim(design)[2] - 1

## Fit the linear model
vfit <-lmFit(corrected_expression, weights = weights, design)
vfit <- eBayes(vfit)

## Create columns for results file
betas = as.data.frame(vfit$coefficients[, ncol(vfit)]); colnames(betas)[1] <- "betas"
p_values = as.data.frame(vfit$p.value[, ncol(vfit)]); colnames(p_values)[1] <- "pvalues"
t_statistic = as.data.frame(vfit$t[, ncol(vfit)]); colnames(t_statistic)[1] <- "t-statistic"
## Uses BH FDR method
fdrs = as.data.frame(p.adjust(p_values[,1], method = "BH", n = length(p_values[,1]))); colnames(fdrs)[1] <- "fdrs"

## Create results df
results <- cbind(betas,t_statistic, p_values, fdrs)


############################################
######## PERMUTE CONDITION LABELS  #########
############################################

# how many interations to do?
iterations = 1000

## Declare permutation function which shuffles the conditions levels.
permute_condition <- function(info, id_colname, trt_colname){
  stopifnot("Treatment column must have exactly two levels."=length(levels(info[,trt_colname]))==2)
  #For each pair of samples
  for(i in 1:length(levels(info[,id_colname]))){
    original_rows <- which(info[,id_colname]==levels(info[,id_colname])[i])
    #If the length of the levels is greater than 1 (in this case it's 2)
    if(length(original_rows)>1){
      #Permute data
      permuted_rows <- sample(original_rows)
      #Write permutation
      info[original_rows,trt_colname] <- info[permuted_rows,trt_colname]
      #Else, shuffle column without sample boundaries
    }else{
      if(runif(1)<0.5){
        # This is NI condition
        info[original_rows,trt_colname]<-levels(info[,trt_colname])[1]
        # This is flu condition
      }else{info[original_rows,trt_colname]<-levels(info[,trt_colname])[2]}
    }
  }
  return(info)
}

###### Do permutations ######
flipped<-c()

if(reproduce){
  ### to reproduce results previously generated load the permuted data.
  shuffled_pvals <- read.table(paste0(results_dir,"/permuted_p_values/permDE.txt"))
  ##if reproduce is set to TRUE this part below will be skipped. If reproduce = FALSE then you the code below will permute the data.
}else{
  meta_data_for_permutations=meta_data
    for(iter in 1:iterations)
    {
      print(iter)

      meta_data_for_permutations_post <-  permute_condition(info=meta_data_for_permutations, id_colname="Genotyping_ID", trt_colname="Condition")
      ## run DE model with permuted conditions
      design = model.matrix(~0 + Genotyping_ID + Condition, data = meta_data_for_permutations_post)
      design <- design[, colSums(design != 0) > 0]
      vfit <-lmFit(corrected_expression, weights = weights, design)
      vfit <- eBayes(vfit)
      permuted_pvals = as.data.frame(vfit$p.value[, ncol(vfit)]); colnames(p_values)[1] <- "pvalues"
      colnames(permuted_pvals) <- iter
      ## Use to check to make sure permutations worked (Condition is column 3)
      ## Does everything but the condition column match?
      identical(meta_data_for_permutations[, names(meta_data_for_permutations) != "Condition"],meta_data_for_permutations_post[, names(meta_data_for_permutations) != "Condition"])
      ## Does the condition column match?
      identical(meta_data_for_permutations[,"Condition"],meta_data_for_permutations_post[,"Condition"])
      

      ## save pvalues from the permutations
      if(iter==1)
      {
        shuffled_pvals <-data.frame(x=permuted_pvals)
        rownames(shuffled_pvals)=rownames(permuted_pvals)
      } else {
        shuffled_pvals <- cbind(shuffled_pvals,x=permuted_pvals)
      }

      ## See what percentage of permutations flipped conditions (ideally around 50% of samples).
      comp<- c()
      for(i in 1:length(meta_data_for_permutations[,"Condition"])){
        if (i==1)
        {
          comp <-identical(meta_data_for_permutations[i,"Condition"],meta_data_for_permutations_post[i,"Condition"])
        } else{
          comp <- cbind(comp, x=identical(meta_data_for_permutations[i,"Condition"],meta_data_for_permutations_post[i,"Condition"]))
        }
      }

      if(iter==1)
      {
        flipped <- table(comp)["FALSE"]
      } else {
        flipped <- cbind(flipped, table(comp)["FALSE"])
        }
    }

  ## Write number of times condition is flipped.
  ## We want this to be ~50%. Should never be equal to 0 because that means you are modeling the real data. 
  write.table(flipped/dim(meta_data)[1],paste0(results_dir,"/permuted_p_values/flipped_counts.txt"))

  ## Save histogram of flipped distribution.
  pdf(paste0(results_dir,"/permuted_p_values/flipped_counts_histogram.pdf"))
  hist(flipped/dim(meta_data)[1], breaks=5, main="Histogram of shuffled condition labels", xlab="shuffled labels")
  dev.off()

  ## Save flipped distribution range and total samples.
  sink(paste0(results_dir,"/permuted_p_values/flipped_counts_stats.txt"))
  cat(paste0("avg=", mean(flipped)/dim(meta_data)[1], "\n", "min=", min(flipped), "\n", "max=", max(flipped), "\n", "Total samples=", nrow(meta_data_for_permutations)))
  sink()
}

## Create directory for permutations
system(paste0("mkdir -p " , results_dir,"/permuted_p_values"))
## Write permutation tests p-values tables
write.table(shuffled_pvals,paste0(results_dir,"/permuted_p_values/permDE.txt"))


###########################################
## SUMMARIZE AND PERFORM FDR CORRECTION  ##
###########################################

## Calculate qvalue
emp_pvals <- empPvals(stat=-log10(results$pvalues), stat0=-log10(as.matrix(shuffled_pvals)), pool = T)
qvals<- qvalue(p=emp_pvals)$qvalue

## Bind qvalue to results tables.
resultsALL=cbind(results, qvals)

## Write results file
assign(paste0("results_DE_infection"), results)
write.table(resultsALL, paste0(results_dir,"/resultsALL.txt"), quote = FALSE, sep = ",")

numGenes_total <- dim(dge)[1]
numDEgenes_10 <- length(which(fdrs[,1] < 0.10))
numDEgenes_05 <- length(which(fdrs[,1] < 0.05))
numDEgenes_01 <- length(which(fdrs[,1] < 0.01))
numDEgenes_001 <- length(which(fdrs[,1] < 0.001))
numDEgenes_0001 <- length(which(fdrs[,1] < 0.0001))

## Create df with results
BHsummary <- cbind(numGenes_total, numDEgenes_10, numDEgenes_05, numDEgenes_01,numDEgenes_001,numDEgenes_0001)

numDEgenes_q10 <- length(which(resultsALL$qvals < 0.10))
numDEgenes_q05 <- length(which(resultsALL$qvals  < 0.05))
numDEgenes_q01 <- length(which(resultsALL$qvals  < 0.01))
numDEgenes_q001 <- length(which(resultsALL$qvals  < 0.001))
numDEgenes_q0001 <- length(which(resultsALL$qvals  < 0.0001))

## Create df with results
qvalsummary<- cbind(numGenes_total, numDEgenes_q10, numDEgenes_q05, numDEgenes_q01, numDEgenes_q001, numDEgenes_q0001)

## Bind these two together
summary <- rbind(BHsummary, qvalsummary)
rownames(summary) <- c("BHsummary", "qvalsummary")

## Write out total number of genes for reference.
write.table(summary, paste0(results_dir,"/summary_info.txt"), quote = FALSE, sep = ",")

## Plot histograms of pvalues, qvalues and B-H FDR and permuted pvalues for QC.
pdf(paste0(results_dir,"/histograms.pdf"))

par(mfrow=c(2,2))
hist(results$pvalues, main="Histogram of raw pvalues")
hist(as.matrix(shuffled_pvals), main="Permuted pvalues")
hist(qvals)
hist(fdrs[,1], main="Histogram of B-H FDR (no perm)")
dev.off()
