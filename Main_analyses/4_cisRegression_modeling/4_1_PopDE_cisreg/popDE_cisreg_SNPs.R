## Load libraries
library(limma)
library(edgeR)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(mgsub)
library(qvalue)

#################################
##### LOAD COMMAND LINE ARGS ####
#################################

args<-commandArgs(TRUE)
#datatype (for directory structure)
data <- args[1]

if(length(args)==0)
  {
    print("WARNING: No arguments supplied.")
  }
print(args)

# Perform 5 permutations to use for FDR correction.
iterations <- 5
reproduce <- "FALSE"


##############################################
## Create directory structure and load data ##
##############################################

folder = "4_cisregression_modeling"
## Set directory 3 steps above script.
setwd('../../../')

## Create directory structure to save outputs.
system(paste0("mkdir -p Outputs/",folder,"/cisregSNPs/", data,"/"))
out_dir <- paste0("Outputs/",folder,"/cisregSNPs/", data,"/")

#Set inputs

##Load metadata
cols = read.table(paste0("Inputs/metadata/", data, "_metadata.txt"))
## Load age and batch corrected voom normalized read counts
reads_whole = read.table(paste0("Outputs/2_ancestry_effects_modeling/batch_corrected_cts_matrices/",data, "_batch.age.corrected.txt"), header=TRUE, row.names=1)

# remove X, Y, MT and contig features from the matrix.
positions <- read.table(paste0("Inputs/QTL_mapping/", data, "_positions.txt"))
autosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
positions_filtered <- positions[positions$chromosome %in% autosomes,]
reads <- reads_whole[rownames(reads_whole) %in% positions_filtered$Gene_ID, ]

## read in genotypes (used in matrixeqtl)
genotypes = read.table(paste0("Inputs/QTL_mapping/SNP_genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)
genotypes$snp <- rownames(genotypes)

## FOR NI condition
## read in top SNPs
NI_topSNPs <- read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/NI/", data,"_NI/results_best_SNPs_with_qval.txt"))
## only keep snp and gene columns
NI_topSNPs <- select(NI_topSNPs, snps, gene)
## subset genotypes to pertinent ones
NI_genotypes_subset <- right_join(genotypes, NI_topSNPs, by = c("snp" = "snps"))
rownames(NI_genotypes_subset) <- NI_genotypes_subset$gene; NI_genotypes_subset$gene <- NULL; NI_genotypes_subset$snp <- NULL
## Change missing data from -9 to NA.
NI_genotypes_subset[NI_genotypes_subset == -9] <- NA

## prepare for adding genotype to meta data
NI_gene_genotype <- as.data.frame(t(NI_genotypes_subset[rownames(NI_genotypes_subset),]))
## add NI to row name
rownames(NI_gene_genotype) <- paste0(rownames(NI_gene_genotype), "_NI")

## FOR FLU condition
## read in top SNPs
FLU_topSNPs <- read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/Flu/", data,"_Flu/results_best_SNPs_with_qval.txt"))
## only keep snp and gene columns
FLU_topSNPs <- select(FLU_topSNPs, snps, gene)
## subset genotypes to pertinent ones
FLU_genotypes_subset <- right_join(genotypes, FLU_topSNPs,  by = c("snp" = "snps"))
rownames(FLU_genotypes_subset) <- FLU_genotypes_subset$gene; FLU_genotypes_subset$gene <- NULL; FLU_genotypes_subset$snp <- NULL
## Change missing data from -9 to NA.
FLU_genotypes_subset[FLU_genotypes_subset == -9] <- NA

## prepare for adding genotype to meta data
FLU_gene_genotype <- as.data.frame(t(FLU_genotypes_subset[rownames(FLU_genotypes_subset),]))
## add NI to row name
rownames(FLU_gene_genotype) <- paste0(rownames(FLU_gene_genotype), "_Flu")

## merge
NI_gene_genotype <- t(NI_gene_genotype)
NI_gene_genotype <- as.data.frame(NI_gene_genotype)

FLU_gene_genotype <- t(FLU_gene_genotype)
FLU_gene_genotype <- as.data.frame(FLU_gene_genotype)

NI_gene_genotype$feature <- rownames(NI_gene_genotype)
FLU_gene_genotype$feature <- rownames(FLU_gene_genotype)

both_gene_genotype <-full_join(NI_gene_genotype, FLU_gene_genotype, by="feature")
rownames(both_gene_genotype) <- both_gene_genotype$feature; both_gene_genotype$feature <- NULL


######################################
## Format reads matrices & metadata ##
######################################

#subset to only include those samples for which there is a 1 in the "PopDEset" column (gives you an opportunity to exclude any samples). (need to add this to metadata files.)
cols=cols[which(cols$PopDE_set==1),]
#ensure that the Sample_ID are the rownames of the cols dataframe
rownames(cols)=cols$Sample
#ensure that the cols dataframe is ordered alphabetically (based on sample_id)
cols=cols[order(rownames(cols)),]

rownames(cols)<- sub("_H3K27ac", "", rownames(cols))
rownames(cols)<- sub("_H3K4me1", "", rownames(cols))
rownames(cols)<- sub("_H3K27me3", "", rownames(cols))
rownames(cols)<- sub("_H3K4me3", "", rownames(cols))

#make sure that all the columns of the meta data are also in the read counts matrix.
reads=reads[,which(colnames(reads) %in% rownames(cols))]
#order the reads count matrix.
reads=reads[,order(colnames(reads))]

#this should be 0. It is a check to ensure that all the metadata and read counts matrix have all the same samples.
length(which(rownames(cols)!=colnames(reads)))
# 0

################################################
## Model DE within each condition (real data) ##
################################################

## subset on features that were originally popDE
## load those popde results
orig_popDE_NI <- read.table(paste0("Outputs/2_ancestry_effects_modeling/PopDE_analysis/", data, "/results/popDE_NI.txt"))
orig_popDE_FLU <- read.table(paste0("Outputs/2_ancestry_effects_modeling/PopDE_analysis/", data, "/results/popDE_FLU.txt"))

## get sig.
sig_NI <- orig_popDE_NI[orig_popDE_NI$qvals_NI < .10, ]
sig_Flu <- orig_popDE_FLU[orig_popDE_FLU$qvals_Flu < .10, ]

## get features that are sig. in either
sig_in_either <- c(rownames(sig_NI), rownames(sig_Flu))
sig_in_either <- unique(sig_in_either)
length(sig_in_either)

## filter down read counts matrix to those that are sig in either
reads <- reads[rownames(reads) %in% sig_in_either, ]

##################################
#### Add genotype to metadata  ###
##################################

sig_genes <- FLU_genotypes_subset[rownames(FLU_genotypes_subset) %in% sig_in_either, ]
genes <- rownames(sig_genes)

for(j in 1:length(genes)){
  if(j%%10 == 0) print(j)
  gene <- genes[j]

    ## add genotype into meta data
  gene_genotype <- as.data.frame(t(both_gene_genotype[rownames(both_gene_genotype) %in% gene,]))
  colnames(gene_genotype)[1] <- "genotype"
  gene_genotype$indiv_ID <- rownames(gene_genotype)
  ## should be 0.
  length(which(rownames(gene_genotype) != gene_genotype$indiv_ID))
  meta_data_loop <- merge(cols, gene_genotype, by=0)
  meta_data_loop$genotype <- as.numeric(meta_data_loop$genotype)
  rownames(meta_data_loop) <- meta_data_loop$Row.names;  meta_data_loop$Row.names <- NULL

  ########################
  ## WITH NESTED MODEL  ##
  ########################

  DE_WC=function(reads,cols){
    design = model.matrix(~0+ Condition + Condition:Admixture + Condition:genotype, data = cols)
    # Remove mock columns
    design <- design[, colSums(design != 0) > 0]
    # Remove individuals with an NA genotype
    reads_loop <- reads[,colnames(reads) %in% rownames(design)]
    fit <-lmFit(reads_loop,design)
    fit <- eBayes(fit)
    return(list(fit,reads_loop))
  }

  DE_BOTH=DE_WC(reads=reads,cols=meta_data_loop)
  fit_BOTH=DE_BOTH[[1]]
  v_BOTH=DE_BOTH[[2]]

  ## for each gene, store each output separately
  results_NI=topTable(fit_BOTH[rownames(fit_BOTH) %in% gene , ],coef="ConditionNI:Admixture",number=nrow(fit_BOTH$coefficients))[,1:4]
  results_FLU=topTable(fit_BOTH[rownames(fit_BOTH) %in% gene , ],coef="ConditionFlu:Admixture",number=nrow(fit_BOTH$coefficients))[,1:4]

  if(j == 1){
    NI_outs <- results_NI
    FLU_outs <- results_FLU
  }else{
    NI_outs <- rbind(NI_outs, results_NI)
    FLU_outs <- rbind(FLU_outs, results_FLU)
  }
}

results_NI <- NI_outs
results_FLU <- FLU_outs

write.table(results_NI, paste0(out_dir, "popDE_NI.real.txt"))
write.table(results_FLU, paste0(out_dir, "popDE_FLU.real.txt"))


#########################################################
## Do permutation tests to define null model for PopDE ##
#########################################################

if(reproduce)
{
  ### to reproduce results previously generated load the permuted data.
  shuffled_pvals_adm_NI=read.table(paste0(out_dir, "/permuted_pvalues/popDE_NI.txt"))
  shuffled_pvals_adm_FLU=read.table(paste0(out_dir, "/permuted_pvalues/popDE_FLU.txt"))
  ##if reproduce is set to TRUE this part below will be skipped. If reproduce = FALSE then you the code below will permute the data.

}else{

  cols_random=cols

  for(iter in 1:iterations)
  {
    if(iter%%100==0)print(iter)
    cols_random$Admixture=sample(cols_random$Admixture)

    for(j in 1:length(genes)){
      if(j%%1000 == 0) print(j)
      gene <- genes[j]

      ## add genotype into meta data
      gene_genotype <- as.data.frame(t(both_gene_genotype[rownames(both_gene_genotype) %in% gene,]))
      colnames(gene_genotype)[1] <- "genotype"
      gene_genotype$indiv_ID <- rownames(gene_genotype)
      ## should be 0.
      length(which(rownames(gene_genotype) != gene_genotype$indiv_ID))
      meta_data_loop <- merge(cols_random, gene_genotype, by=0)
      meta_data_loop$genotype <- as.numeric(meta_data_loop$genotype)
      rownames(meta_data_loop) <- meta_data_loop$Row.names;  meta_data_loop$Row.names <- NULL

      ########################
      ## WITH NESTED MODEL  ##
      ########################

      DE_WC=function(reads,cols){
        design = model.matrix(~0+ Condition + Condition:Admixture + Condition:genotype, data = cols)
        # Remove mock columns
        design <- design[, colSums(design != 0) > 0]
        # Remove individuals with an NA genotype
        reads_loop <- reads[,colnames(reads) %in% rownames(design)]
        fit <-lmFit(reads_loop,design)
        fit <- eBayes(fit)
        return(list(fit,reads_loop))
      }

      DE_BOTH=DE_WC(reads=reads,cols=meta_data_loop)
      fit_BOTH=DE_BOTH[[1]]
      v_BOTH=DE_BOTH[[2]]

      ## for each gene, store each output separately
      perm_NI=topTable(fit_BOTH[rownames(fit_BOTH) %in% gene , ],coef="ConditionNI:Admixture",number=nrow(fit_BOTH$coefficients))[,1:4]
      perm_FLU=topTable(fit_BOTH[rownames(fit_BOTH) %in% gene , ],coef="ConditionFlu:Admixture",number=nrow(fit_BOTH$coefficients))[,1:4]

      if(j == 1){
        NI_outs <- perm_NI
        FLU_outs <- perm_FLU
      }else{
        NI_outs <- rbind(NI_outs, perm_NI)
        FLU_outs <- rbind(FLU_outs, perm_FLU)
      }
    }
    
    permuted_adm_FLU=FLU_outs
    permuted_adm_NI=NI_outs

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
  write.table(shuffled_pvals_adm_FLU,paste0(out_dir, "/permuted_pvalues/popDE_FLU.txt"))
  write.table(shuffled_pvals_adm_NI,paste0(out_dir, "/permuted_pvalues/popDE_NI.txt"))

}


########################################################################
## Write results and add additional FDR correction methods to results ##
########################################################################

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
write.table(results_FLU_BH_q, paste0(out_dir, "popDE_FLU.qval.txt"))
write.table(results_NI_BH_q, paste0(out_dir, "popDE_NI.qval.txt"))

############################
## Print sig. PopDE genes ##
############################

#write this as an output and place in summary directory.
sink(paste0(out_dir, "summary_info.txt"))

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