## Load libraries
library(dplyr)
library(relaimpo)

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
reads_whole = read.table(paste0("Outputs/1_DE_infection_modeling/", data, "_residuals/corrected_expression.txt"), header = TRUE, sep = ",", row.names = 1)

# remove X, Y, MT and contig features from the matrix.
positions <- read.table(paste0("Inputs/QTL_mapping/", data, "_positions.txt"))
autosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
               "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
positions_filtered <- positions[positions$chromosome %in% autosomes,]
reads <- reads_whole[rownames(reads_whole) %in% positions_filtered$Gene_ID, ]


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

## filter down read counts matrix to those that are sig in either
#reads <- reads[rownames(reads) %in% sig_in_either, ]
reads <- as.matrix(reads)

##################################
#### Add genotype to metadata  ###
##################################

both_variants_genes <- intersect(rownames(FLU_STRgenotypes_subset), rownames(FLU_SNPgenotypes_subset))
#genes <- both_variants_genes[which(both_variants_genes %in% sig_in_either)]
genes <- both_variants_genes

results_FLU = data.frame(matrix(nrow = length(genes), ncol = 4))
rownames(results_FLU) <- genes
colnames(results_FLU) <- c("no_variant", "SNP", "STR", "SNP_STR")

results_NI = data.frame(matrix(nrow = length(genes), ncol = 4))
rownames(results_NI) <- genes
colnames(results_NI) <- c("no_variant", "SNP", "STR", "SNP_STR")

for(j in 1:length(genes)){
  if(j%%10 == 0) print(j)
  gene <- genes[j]
  
  meta_data_loop <- cols
  meta_data_loop$indiv_ID <- rownames(meta_data_loop)
  
  ## add genotype into meta data
  meta_data_loop$SNP <- as.numeric(both_gene_SNPgenotype[gene, rownames(meta_data_loop)])
  meta_data_loop$STR <- as.numeric(both_gene_STRgenotype[gene, rownames(meta_data_loop)])
  
  Flu_index <-seq(1,ncol(reads),2)
  reads_Flu <- reads[,Flu_index]
  NI_index <-seq(2,ncol(reads),2)
  reads_NI <- reads[,NI_index]
  
  meta_data_loop_NI <- subset(meta_data_loop, meta_data_loop$Condition=="NI")
  meta_data_loop_Flu <- subset(meta_data_loop, meta_data_loop$Condition=="Flu")
  
  ###########################
  ## modeling and relaimpo ##
  ###########################
  
  # model without variant genotype
  gene_model_Flu=lm(reads_Flu[gene,] ~ Age + Admixture, data = meta_data_loop_Flu)
  tryCatch(
    {
      rel_impo_Flu=suppressWarnings(calc.relimp(gene_model_Flu,type=c("last")))
      results_FLU[j,1]=rel_impo_Flu$last["Admixture"]  
    }, 
    error=function(cond) {
      message(cond)
      message(" - SKIPPED")
      results_FLU[j,1]="NA"
    }
  )  
  gene_model_NI=lm(reads_NI[gene,] ~ Age + Admixture, data = meta_data_loop_NI)
  tryCatch(
    {
      rel_impo_NI=suppressWarnings(calc.relimp(gene_model_NI,type=c("last")))
      results_NI[j,1]=rel_impo_NI$last["Admixture"]
    }, 
    error=function(cond) {
      message(cond)
      message(" - SKIPPED")
      results_NI[j,1]="NA"
    }
  )   
  
  # model with SNP genotype
  gene_model_Flu=lm(reads_Flu[gene,] ~ Age + Admixture + SNP, data = meta_data_loop_Flu)
  tryCatch(
    {
      rel_impo_Flu=suppressWarnings(calc.relimp(gene_model_Flu,type=c("last")))
      results_FLU[j,2]=rel_impo_Flu$last["Admixture"]
    }, 
    error=function(cond) {
      message(cond)
      message(" - SKIPPED")
      results_FLU[j,2]="NA"
    }
  )  
  gene_model_NI=lm(reads_NI[gene,] ~ Age + Admixture + SNP, data = meta_data_loop_NI)
  tryCatch(
    {
      rel_impo_NI=suppressWarnings(calc.relimp(gene_model_NI,type=c("last")))    
      results_NI[j,2]=rel_impo_NI$last["Admixture"]
    }, 
    error=function(cond) {
      message(cond)
      message(" - SKIPPED")
      results_NI[j,2]="NA"
    }
  )  
    
  # model with STR genotype
  gene_model_Flu=lm(reads_Flu[gene,] ~ Age + Admixture + STR, data = meta_data_loop_Flu)
  tryCatch(
    {
      rel_impo_Flu=suppressWarnings(calc.relimp(gene_model_Flu,type=c("last")))
      results_FLU[j,3]=rel_impo_Flu$last["Admixture"]
    }, 
    error=function(cond) {
      message(cond)
      message(" - SKIPPED")
      results_FLU[j,3]="NA"
    }
  )  
  
  gene_model_NI=lm(reads_NI[gene,] ~ Age + Admixture + STR, data = meta_data_loop_NI)
  tryCatch(
    {
      rel_impo_NI=suppressWarnings(calc.relimp(gene_model_NI,type=c("last")))    
      results_NI[j,3]=rel_impo_NI$last["Admixture"]  
    }, 
    error=function(cond) {
      message(cond)
      message(" - SKIPPED")
      results_NI[j,3]="NA"
    }
  )    
  
  # model with SNP and STR genotypes
  gene_model_Flu=lm(reads_Flu[gene,] ~ Age + Admixture + SNP + STR, data = meta_data_loop_Flu)
  tryCatch(
    {
      rel_impo_Flu=suppressWarnings(calc.relimp(gene_model_Flu,type=c("last")))
      results_FLU[j,4]=rel_impo_Flu$last["Admixture"]  
    }, 
    error=function(cond) {
      message(cond)
      message(" - SKIPPED")
      results_FLU[j,4]="NA"
    }
  )
  
  gene_model_NI=lm(reads_NI[gene,] ~ Age + Admixture + SNP + STR, data = meta_data_loop_NI)
  tryCatch(
    {
      rel_impo_NI=suppressWarnings(calc.relimp(gene_model_NI,type=c("last")))    
      results_NI[j,4]=rel_impo_NI$last["Admixture"]
    }, 
    error=function(cond) {
      message(cond)
      message(" - SKIPPED")
      results_NI[j,4]="NA"
    }
  )      
}



write.table(results_NI, paste0(out_dir, "all_feature_rel_impo_NI_last.txt"))
write.table(results_FLU, paste0(out_dir, "all_feature_rel_impo_FLU_last.txt"))
