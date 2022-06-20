

library(bsseq)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(qvalue)
library(DSS)

library(relaimpo)

#################################
##### LOAD COMMAND LINE ARGS ####
#################################


#datatype (for directory structure)
data <- "WGBS"


##############################################
## Create directory structure and load data ##
##############################################

folder = "4_cisregression_modeling"
## Set directory 3 steps above script.
setwd('../../../')

## Create directory structure to save outputs.
system(paste0("mkdir -p Outputs/",folder,"/cisregSNPsSTRs/", data,"/"))
out_dir <- paste0("Outputs/",folder,"/cisregSNPsSTRs/", data,"/")

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
#Set inputs
#Load metadata
cols = read.table(paste0("Inputs/metadata/", data, "_metadata.txt"))

## Load unsmoothed BSseq object
load("Inputs/counts_matrices/WGBS_filtered.counts.BSobj.RData")

## Load and subset SNP genotypes
## read in SNP genotypes (used in matrixeqtl)
SNPgenotypes = read.table(paste0("Inputs/QTL_mapping/SNP_genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)
SNPgenotypes$snp <- rownames(SNPgenotypes)

## read in top SNPs for NI condition
NI_topSNPs <- read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/NI/", data,"_NI/results_best_SNPs_with_qval.txt"))
## only keep snp and gene columns
NI_topSNPs <- NI_topSNPs[,c("snps", "gene")]
## subset genotypes to pertinent ones
NI_SNPgenotypes_subset <- right_join(SNPgenotypes, NI_topSNPs, by = c("snp" = "snps"))
rownames(NI_SNPgenotypes_subset) <- NI_SNPgenotypes_subset$gene; NI_SNPgenotypes_subset$gene <- NULL; NI_SNPgenotypes_subset$snp <- NULL
## Change missing data from -9 to NA.
NI_SNPgenotypes_subset[NI_SNPgenotypes_subset == -9] <- NA

## prepare for adding genotype to meta data
NI_gene_SNPgenotype <- NI_SNPgenotypes_subset[,colnames(NI_SNPgenotypes_subset)]
## add NI to col name
colnames(NI_gene_SNPgenotype) <- paste0(colnames(NI_gene_SNPgenotype), "_NI")

## read in top SNPs for FLU condition
FLU_topSNPs <- read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/Flu/", data,"_Flu/results_best_SNPs_with_qval.txt"))
## only keep snp and gene columns
FLU_topSNPs <- FLU_topSNPs[, c("snps", "gene")]
## subset genotypes to pertinent ones
FLU_SNPgenotypes_subset <- right_join(SNPgenotypes, FLU_topSNPs,  by = c("snp" = "snps"))
rownames(FLU_SNPgenotypes_subset) <- FLU_SNPgenotypes_subset$gene; FLU_SNPgenotypes_subset$gene <- NULL; FLU_SNPgenotypes_subset$snp <- NULL
## Change missing data from -9 to NA.
FLU_SNPgenotypes_subset[FLU_SNPgenotypes_subset == -9] <- NA

## prepare for adding genotype to meta data
FLU_gene_SNPgenotype <- FLU_SNPgenotypes_subset[,colnames(FLU_SNPgenotypes_subset)]
## add FLU to row name
colnames(FLU_gene_SNPgenotype) <- paste0(colnames(FLU_gene_SNPgenotype), "_Flu")

## merge the genotypes of the best SNP in NI and Flu condition  
NI_gene_SNPgenotype$feature <- rownames(NI_gene_SNPgenotype)
FLU_gene_SNPgenotype$feature <- rownames(FLU_gene_SNPgenotype)

both_gene_SNPgenotype <-full_join(NI_gene_SNPgenotype, FLU_gene_SNPgenotype, by="feature")
rownames(both_gene_SNPgenotype) <- both_gene_SNPgenotype$feature; both_gene_SNPgenotype$feature <- NULL


## Load and subset STR genotypes
## read in STR genotypes (used in matrixeqtl)
STRgenotypes = read.table(paste0("Inputs/QTL_mapping/STR_genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)
STRgenotypes$str <- rownames(STRgenotypes)

## read in top STRs for NI condition
NI_topSTRs <- read.table(paste0("Outputs/3_QTL_mapping/STR-QTL_mapping/",data,"/NI/", data,"_NI/results_best_STRs_with_qval.txt"))
## only keep str and gene columns
NI_topSTRs <- NI_topSTRs[, c("snps", "gene")]
## subset genotypes to pertinent ones
NI_STRgenotypes_subset <- right_join(STRgenotypes, NI_topSTRs, by = c("str" = "snps"))
rownames(NI_STRgenotypes_subset) <- NI_STRgenotypes_subset$gene; NI_STRgenotypes_subset$gene <- NULL; NI_STRgenotypes_subset$str <- NULL
## Change missing data from -9 to NA.
NI_STRgenotypes_subset[NI_STRgenotypes_subset == -9] <- NA

## prepare for adding genotype to meta data
NI_gene_STRgenotype <- NI_STRgenotypes_subset[,colnames(NI_STRgenotypes_subset)]
## add NI to col name
colnames(NI_gene_STRgenotype) <- paste0(colnames(NI_gene_STRgenotype), "_NI")

## read in top STRs for FLU condition
FLU_topSTRs <- read.table(paste0("Outputs/3_QTL_mapping/STR-QTL_mapping/",data,"/Flu/", data,"_Flu/results_best_STRs_with_qval.txt"))
## only keep snp and gene columns
FLU_topSTRs <- FLU_topSTRs[, c("snps", "gene")]
## subset genotypes to pertinent ones
FLU_STRgenotypes_subset <- right_join(STRgenotypes, FLU_topSTRs,  by = c("str" = "snps"))
rownames(FLU_STRgenotypes_subset) <- FLU_STRgenotypes_subset$gene; FLU_STRgenotypes_subset$gene <- NULL; FLU_STRgenotypes_subset$str <- NULL
## Change missing data from -9 to NA.
FLU_STRgenotypes_subset[FLU_STRgenotypes_subset == -9] <- NA

## prepare for adding genotype to meta data
FLU_gene_STRgenotype <- FLU_STRgenotypes_subset[,colnames(FLU_STRgenotypes_subset)]
## add FLU to col name
colnames(FLU_gene_STRgenotype) <- paste0(colnames(FLU_gene_STRgenotype), "_Flu")

## merge the genotypes of the best STR in NI and Flu condition 
NI_gene_STRgenotype$feature <- rownames(NI_gene_STRgenotype)
FLU_gene_STRgenotype$feature <- rownames(FLU_gene_STRgenotype)

both_gene_STRgenotype <-full_join(NI_gene_STRgenotype, FLU_gene_STRgenotype, by="feature")
rownames(both_gene_STRgenotype) <- both_gene_STRgenotype$feature; both_gene_STRgenotype$feature <- NULL


######################################
## Format reads matrices & metadata ##
######################################

#subset to only include those samples for which there is a 1 in the "PopDEset" column
cols=cols[which(cols$PopDE_set==1),]
#ensure that the Sample_ID are the rownames of the cols dataframe
rownames(cols)=cols$Sample
#ensure that the cols dataframe is ordered alphabetically (based on sample_id)
cols=cols[order(rownames(cols)),]
## Mean center age & refactorize
cols=pretty_up_cols(cols)
## Set NI as first level so that it is used as the reference.
cols$Condition = factor(cols$Condition, levels=c("NI","Flu"))

## Check to make sure the number of samples in metadata and WGBS reads data matches.
dim(cols)[1] == dim(BSobj.fit.nolow)[2]
dim(BSobj.fit.nolow)[1]==7463164 # number of cpg sites there should be


## Filter based on coverage
#To minimize noise in methylation estimates due to low-coverage data, we restricted analysis to CpG sites
# with coverage of ≥4 sequence reads in at least half of the samples in each condition. DSS accounts for low coverage
# in the model so it is not necessary before.
BS.cov <- getCoverage(BSobj.fit.nolow)
keepLoci.ex <- which(rowSums(BS.cov[, BSobj.fit.nolow$Condition == "Flu"] >= 4) > 17 &
                       rowSums(BS.cov[, BSobj.fit.nolow$Condition == "NI"] >= 4) > 17 )
length(keepLoci.ex)
BSobj.fit.nolow <- BSobj.fit.nolow[keepLoci.ex,]

## remove sites with no variation
reads_whole <-  getMeth(BSobj.fit.nolow, type="raw")
no.var <- rownames(reads_whole[rowVars(reads_whole, na.rm=TRUE) == 0, ])
reads_whole.lim <- subset(reads_whole, !(rownames(reads_whole) %in% no.var))
#BSobj.fit.nolow.var <- BSobj.fit.nolow[rownames(reads_whole.lim),]

reads <- reads_whole.lim[,which(colnames(reads_whole.lim) %in% rownames(cols))]

# remove X, Y, MT and contig features from the matrix.
positions <- read.table(paste0("Inputs/QTL_mapping/", data, "_positions.txt"))
autosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
positions_filtered <- positions[positions$chromosome %in% autosomes,]
#BS.fit <- BSobj.fit.nolow.var[rownames(BSobj.fit.nolow.var) %in% positions_filtered$Gene_ID, ]
reads <- reads[rownames(reads) %in% positions_filtered$Gene_ID, ]

#################################
## Subset genes to be analyzed ##
#################################

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

## filter down bsseq object that are sig in either
reads <- reads[rownames(reads) %in% sig_in_either, ]
reads <- as.matrix(reads)

save.image(file=paste0(out_dir,"/SNP_STR_admixturePVE_WGBS.RData"))
