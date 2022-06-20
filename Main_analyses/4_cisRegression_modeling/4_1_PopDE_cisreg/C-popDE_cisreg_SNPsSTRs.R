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
system(paste0("mkdir -p Outputs/",folder,"/cisregSNPsSTRs/", data,"/"))
out_dir <- paste0("Outputs/",folder,"/cisregSNPsSTRs/", data,"/")

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

## read in STR genotypes (used in matrixeqtl)
STRgenotypes = read.table(paste0("Inputs/QTL_mapping/STR_genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)
STRgenotypes$str <- rownames(STRgenotypes)

## read in top STRs
## for NI condition
NI_topSTRs <- read.table(paste0("Outputs/3_QTL_mapping/STR-QTL_mapping/",data,"/NI/", data,"_NI/results_best_STRs_with_qval.txt"))
## only keep str and gene columns
NI_topSTRs <- select(NI_topSTRs, snps, gene)
## subset genotypes to pertinent ones
NI_STRgenotypes_subset <- right_join(STRgenotypes, NI_topSTRs, by = c("str" = "snps"))
rownames(NI_STRgenotypes_subset) <- NI_STRgenotypes_subset$gene; NI_STRgenotypes_subset$gene <- NULL; NI_STRgenotypes_subset$str <- NULL
## Change missing data from -9 to NA.
NI_STRgenotypes_subset[NI_STRgenotypes_subset == -9] <- NA

## prepare for adding genotype to meta data
NI_gene_STRgenotype <- NI_STRgenotypes_subset[,colnames(NI_STRgenotypes_subset)]
## add NI to col name
colnames(NI_gene_STRgenotype) <- paste0(colnames(NI_gene_STRgenotype), "_NI")

## read in top STRs
## for FLU condition
FLU_topSTRs <- read.table(paste0("Outputs/3_QTL_mapping/STR-QTL_mapping/",data,"/Flu/", data,"_Flu/results_best_STRs_with_qval.txt"))
## only keep snp and gene columns
FLU_topSTRs <- select(FLU_topSTRs, snps, gene)
## subset genotypes to pertinent ones
FLU_STRgenotypes_subset <- right_join(STRgenotypes, FLU_topSTRs,  by = c("str" = "snps"))
rownames(FLU_STRgenotypes_subset) <- FLU_STRgenotypes_subset$gene; FLU_STRgenotypes_subset$gene <- NULL; FLU_STRgenotypes_subset$str <- NULL
## Change missing data from -9 to NA.
FLU_STRgenotypes_subset[FLU_STRgenotypes_subset == -9] <- NA

## prepare for adding genotype to meta data
FLU_gene_STRgenotype <- FLU_STRgenotypes_subset[,colnames(FLU_STRgenotypes_subset)]
## add FLU to col name
colnames(FLU_gene_STRgenotype) <- paste0(colnames(FLU_gene_STRgenotype), "_Flu")

## merge

NI_gene_STRgenotype$feature <- rownames(NI_gene_STRgenotype)
FLU_gene_STRgenotype$feature <- rownames(FLU_gene_STRgenotype)

both_gene_STRgenotype <-full_join(NI_gene_STRgenotype, FLU_gene_STRgenotype, by="feature")
rownames(both_gene_STRgenotype) <- both_gene_STRgenotype$feature; both_gene_STRgenotype$feature <- NULL


## read in SNP genotypes (used in matrixeqtl)
SNPgenotypes = read.table(paste0("Inputs/QTL_mapping/SNP_genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)
SNPgenotypes$snp <- rownames(SNPgenotypes)

## read in top SNPs
## for NI condition
NI_topSNPs <- read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/NI/", data,"_NI/results_best_SNPs_with_qval.txt"))
## only keep snp and gene columns
NI_topSNPs <- select(NI_topSNPs, snps, gene)
## subset genotypes to pertinent ones
NI_SNPgenotypes_subset <- right_join(SNPgenotypes, NI_topSNPs, by = c("snp" = "snps"))
rownames(NI_SNPgenotypes_subset) <- NI_SNPgenotypes_subset$gene; NI_SNPgenotypes_subset$gene <- NULL; NI_SNPgenotypes_subset$snp <- NULL
## Change missing data from -9 to NA.
NI_SNPgenotypes_subset[NI_SNPgenotypes_subset == -9] <- NA

## prepare for adding genotype to meta data
NI_gene_SNPgenotype <- NI_SNPgenotypes_subset[,colnames(NI_SNPgenotypes_subset)]
## add NI to col name
colnames(NI_gene_SNPgenotype) <- paste0(colnames(NI_gene_SNPgenotype), "_NI")

## read in top SNPs
## for FLU condition
FLU_topSNPs <- read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/Flu/", data,"_Flu/results_best_SNPs_with_qval.txt"))
## only keep snp and gene columns
FLU_topSNPs <- select(FLU_topSNPs, snps, gene)
## subset genotypes to pertinent ones
FLU_SNPgenotypes_subset <- right_join(SNPgenotypes, FLU_topSNPs,  by = c("snp" = "snps"))
rownames(FLU_SNPgenotypes_subset) <- FLU_SNPgenotypes_subset$gene; FLU_SNPgenotypes_subset$gene <- NULL; FLU_SNPgenotypes_subset$snp <- NULL
## Change missing data from -9 to NA.
FLU_SNPgenotypes_subset[FLU_SNPgenotypes_subset == -9] <- NA

## prepare for adding genotype to meta data
FLU_gene_SNPgenotype <- FLU_SNPgenotypes_subset[,colnames(FLU_SNPgenotypes_subset)]
## add FLU to row name
colnames(FLU_gene_SNPgenotype) <- paste0(colnames(FLU_gene_SNPgenotype), "_Flu")

## merge

NI_gene_SNPgenotype$feature <- rownames(NI_gene_SNPgenotype)
FLU_gene_SNPgenotype$feature <- rownames(FLU_gene_SNPgenotype)

both_gene_SNPgenotype <-full_join(NI_gene_SNPgenotype, FLU_gene_SNPgenotype, by="feature")
rownames(both_gene_SNPgenotype) <- both_gene_SNPgenotype$feature; both_gene_SNPgenotype$feature <- NULL


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

both_variants_genes <- intersect(rownames(FLU_STRgenotypes_subset), rownames(FLU_SNPgenotypes_subset))
genes <- both_variants_genes[which(both_variants_genes %in% sig_in_either)]

for(j in 1:length(genes)){
  if(j%%10 == 0) print(j)
  gene <- genes[j]

  meta_data_loop <- cols
  meta_data_loop$indiv_ID <- rownames(meta_data_loop)
  
    ## add genotype into meta data
  meta_data_loop$SNP <- as.numeric(both_gene_SNPgenotype[gene, rownames(meta_data_loop)])
  meta_data_loop$STR <- as.numeric(both_gene_STRgenotype[gene, rownames(meta_data_loop)])

  ########################
  ## WITH NESTED MODEL  ##
  ########################

  DE_WC=function(reads,cols){
    design = model.matrix(~0+ Condition + Condition:Admixture + Condition:SNP + Condition:STR, data = cols)
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
      
      meta_data_loop <- cols_random
      meta_data_loop$indiv_ID <- rownames(meta_data_loop)
      
      ## add genotype into meta data
      meta_data_loop$SNP <- as.numeric(both_gene_SNPgenotype[gene, rownames(meta_data_loop)])
      meta_data_loop$STR <- as.numeric(both_gene_STRgenotype[gene, rownames(meta_data_loop)])
      

      ########################
      ## WITH NESTED MODEL  ##
      ########################

      DE_WC=function(reads,cols){
        design = model.matrix(~0+ Condition + Condition:Admixture + Condition:SNP + Condition:STR, data = cols)
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
