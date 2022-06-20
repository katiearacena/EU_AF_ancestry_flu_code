## Perform DE infection analysis for WGBS

## Load libraries
library(DSS)
library(bsseq)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(qvalue)

##  TRUE if want to perform permutations. FALSE if want to use already created permutations.
reproduce <- TRUE
iterations <- 10
data <- "WGBS"

#############################################
## CREATE DIRECTORY STRUCTURE & LOAD DATA  ##
#############################################
## Set directory structure
folder = "1_DE_infection_modeling"
## Set directory 3 steps above script.
setwd('../../../')
## Set results directory
system(paste0("mkdir -p Outputs/",folder, "/WGBS_results"))
results_dir <- paste0("Outputs/",folder, "/WGBS_results/")

## metadata
meta_data = read.table(paste0("Inputs/metadata/", data, "_metadata.txt"))
## Load filtered BSobj containing WGBS cts data.
load("Inputs/counts_matrices/WGBS_filtered.counts.BSobj.RData")
BS.fit <- BSobj.fit.nolow

#############################
###### Filter metadata ######
#############################

## Remove unneeded rows from the metadata file.
## Remove Mock and IPSC samples
meta_data <- meta_data[which(meta_data$PopDE_set=="1"),]

## Factorize certain columns from meta data
## *Batch MUST be a factor so you get multiple estimates!*
meta_data$Batch <- as.factor(meta_data$Batch)
meta_data$Genotyping_ID <- as.factor(meta_data$Genotyping_ID)
## Set NI as first level so that it is used as the reference.
meta_data$Condition = factor(meta_data$Condition, levels=c("NI","Flu"))

## Check to make sure the number of samples in metadata and WGBS reads data matches.
dim(meta_data)[1] == dim(BS.fit)[2]
dim(BS.fit)[1] # number of cpg sites that pass filters


#############################
## Computing t-statistics ###
#############################

## Set design matrix
design <- pData(BS.fit)
X <- model.matrix(~Condition + pair, data = design)
X

## Fit model using the design matrix.
## NO Smoothing here
DMLfit = DMLfit.multiFactor(BS.fit, design=design, formula= ~Condition + pair, smoothing = FALSE)
dmlTest = DMLtest.multiFactor(DMLfit, term="Condition")
rownames(dmlTest) <- granges(BS.fit)

save(dmlTest, file=paste0(results_dir, "results_of_DSS_dml_test.RData"))
write.table(dmlTest, file=paste0(results_dir, "results_dmlTest.txt"))

pdf(paste0(results_dir, "DSS.allchrs.NO.SMOOTHING.pval_hist.pdf"))
hist(dmlTest$pvals, 200, main="all chrs & no smoothing")
dev.off()

pvals <- as.data.frame(dmlTest$pvals)
rownames(pvals) <- rownames(dmlTest)

## Declare permutation function (from Paul)
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

    ## do test
    pData(BS.fit)$pair <- factor(meta_data_for_permutations_post$Genotyping_ID)
	  pData(BS.fit)$Condition <- factor(meta_data_for_permutations_post$Condition)
    design <- pData(BS.fit)
    DMLfit = DMLfit.multiFactor(BS.fit, design=design, formula= ~ Condition + pair, smoothing = F)
	  permuted.tstats =  DMLtest.multiFactor(DMLfit, term="Condition")
	  rownames(permuted.tstats) <- granges(BS.fit)
  	perm.pvals <- as.data.frame(permuted.tstats$pvals)
	  rownames(perm.pvals) <- rownames(permuted.tstats)

    ## Use to check to make sure permutations worked (Condition is column 3)
    ## Does everything but the condition column match? **
    identical(meta_data_for_permutations[, names(meta_data_for_permutations) != "Condition"],meta_data_for_permutations_post[, names(meta_data_for_permutations) != "Condition"])
    ## Does the condition column match?
    identical(meta_data_for_permutations[,"Condition"],meta_data_for_permutations_post[,"Condition"])


    if(iter==1)
    {
      shuffled_pvals <-data.frame(x=perm.pvals)
      rownames(shuffled_pvals)=rownames(perm.pvals)
      colnames(shuffled_pvals)=iter
    } else {
      shuffled_pvals <- cbind(shuffled_pvals,x=perm.pvals)
    }

    ## See if % of permutations flipped (ideally around 50%).
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
  ## We want this to be ~50% and defintely don't want this value to be 0 because that means that 0 are flipped and it's equal to the real data.
  write.table(flipped,paste0(results_dir,"/permuted_p_values/flipped_counts.txt"))

  ## Save histogram of flipped distribution.
  pdf(paste0(results_dir,"/permuted_p_values/flipped_counts_histogram.pdf"))
  hist(flipped, breaks=5, main="Histogram of shuffled condition labels", xlab="shuffled labels")
  dev.off()

  ## Save flipped distribution range and total samples.
  sink(paste0(results_dir,"/permuted_p_values/flipped_counts_stats.txt"))
  cat(paste0("min=", min(flipped), "\n", "max=", max(flipped), "\n", "Total samples=", nrow(meta_data_for_permutations)))
  sink()
}

## Create directory for permutations
system(paste0("mkdir -p " , results_dir,"/permuted_p_values"))
## Write permutation tests p-values tables
colnames(shuffled_pvals) <- c(1:iterations)
write.table(shuffled_pvals,paste0(results_dir,"/permuted_p_values/permDE.txt"))

### Plot histogram of permuted pvalue histograms
shuffled_pvals %>% gather() %>% head()
pdf(paste0(results_dir,"perm_pvals_hist.pdf"))
ggplot(gather(shuffled_pvals), aes(value)) +
  geom_histogram(bins = 100) +
  facet_wrap(~key, scales = 'free_x')
dev.off()

###########################################
## SUMMARIZE AND PERFORM FDR CORRECTION  ##
###########################################

## Calculate qvalue
pvals <- na.omit(pvals)
shuffled_pvals <- na.omit(shuffled_pvals)

emp_pvals <- empPvals(stat=-log10(as.matrix(pvals)), stat0=-log10(as.matrix(shuffled_pvals)), pool = T)
qvals<- qvalue(p=emp_pvals)$qvalue

## Bind qvalue to results tables.
results <- na.omit(dmlTest)
resultsALL=cbind(results, qvals)

## Add average methylation counts and differences to results file
cts <- getMeth(BS.fit, type="raw")
NI_samples <- cts[,grepl("NI", colnames(cts))]
FLU_samples <- cts[,grepl("Flu", colnames(cts))]
NI_avg <- as.data.frame(rowMeans(NI_samples, na.rm=TRUE))
colnames(NI_avg) <- c("NI_mean")
FLU_avg <- as.data.frame(rowMeans(FLU_samples, na.rm=TRUE))
colnames(FLU_avg) <- c("Flu_mean")
group_means <- merge(NI_avg, FLU_avg, by=0, all=T)
rownames(group_means) <- group_means$Row.names; group_means$Row.names=NULL
## now change colon to underscore
rownames(resultsALL) <- str_replace(rownames(resultsALL), ":", "_")

## now merge with dml results
resultsALL$ID <- rownames(resultsALL)
group_means$ID <- rownames(group_means)
results_gmeans <- left_join(resultsALL, group_means, by="ID")
rownames(results_gmeans) <- results_gmeans$ID; results_gmeans$ID <- NULL

results_gmeans$diff <- results_gmeans$NI_mean - results_gmeans$Flu_mean

## Write results file
write.table(results_gmeans, paste0(results_dir, "resultsALL.txt"), quote = FALSE, sep = ",")

numGenes_total <- dim(results_gmeans)[1]
numDEgenes_q20 <- length(which(results_gmeans$qvals < 0.20))
numDEgenes_q10 <- length(which(results_gmeans$qvals < 0.10))
numDEgenes_q20_md <- length(which(results_gmeans$qvals < 0.20 & (results_gmeans$diff >= .10 | results_gmeans$diff <= -.10)))
numDEgenes_q10_md <- length(which(results_gmeans$qvals < 0.10 & (results_gmeans$diff >= .10 | results_gmeans$diff <= -.10)))

## Create df with results
qvalsummary<- cbind(numGenes_total, numDEgenes_q20, numDEgenes_q10, numDEgenes_q20_md, numDEgenes_q10_md)

## Write out total number of genes for reference.
write.table(qvalsummary, paste0(results_dir, "summary_info.txt"), quote = FALSE, sep = ",")