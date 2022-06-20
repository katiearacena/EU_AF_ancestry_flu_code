## Load libraries
library(limma)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(data.table)
library(gridExtra)
library(mgsub)
library(metap)
library(qusage)
library(msigdbr)
library(biomaRt)
library(org.Hs.eg.db)
library(bsseq)
library(ggridges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(tidyverse)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(GSVA)
library(data.table)

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

## Get human pathways
m_df = msigdbr(species = "Homo sapiens")
## Check the available collections and sub-collections
m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
## Retrieve human genes for the hallmark collection gene sets
m_df = msigdbr(species = "Homo sapiens", category = "H")
human.path.list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

## select pathway you want gene list for
gene_list_tfna <- human.path.list$HALLMARK_TNFA_SIGNALING_VIA_NFKB
gene_list_inflam <- human.path.list$HALLMARK_INFLAMMATORY_RESPONSE
gene_list_IL2 <- human.path.list$HALLMARK_IL2_STAT5_SIGNALING
gene_list_IL6 <- human.path.list$HALLMARK_IL6_JAK_STAT3_SIGNALING
gene_list_ifna<- human.path.list$HALLMARK_INTERFERON_ALPHA_RESPONSE
gene_list_ifng <- human.path.list$HALLMARK_INTERFERON_GAMMA_RESPONSE

## get union of 5 inflammatory pathways (excluding ifna)
combo1 <- Reduce(union, list(gene_list_tfna, gene_list_inflam, gene_list_IL2, gene_list_IL6, gene_list_ifng))
combo2 <- Reduce(union, list(gene_list_tfna, gene_list_inflam, gene_list_IL2, gene_list_IL6))

## create list of all gene sets you want to test

gs2 <- list(gene_list_tfna, gene_list_inflam, gene_list_IL2, gene_list_IL6, gene_list_ifna, gene_list_ifng, combo1, combo2)
names(gs2) <- c("TFNA", "INFLAM", "IL2", "IL6", "IFNA", "IFNG", "TFNA_INFLA_IL2_IL6_IFNG", "TFNA_INFLA_IL2_IL6")


pathway_genes <- unique(unlist(gs2)) 
##number of unique immune genes
length(pathway_genes)

#############################################
## CREATE DIRECTORY STRUCTURE & LOAD DATA  ##
#############################################

## Set directory structure
folder = "4_cisregression_modeling"
## Set directory 3 steps above script.
setwd('../../../')

system(paste0("mkdir -p Outputs/",folder,"/cisreg_ancestry_score/", data,"/"))
out_dir <- paste0("Outputs/",folder,"/cisreg_ancestry_score/", data,"/")

## REGRESS OUT THE TOP SNP AND STR FROM COUNTS MATRIX

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


## Get gene-peak pairs info
if(data=="RNAseq"){
	results_with_genes <- read.table(paste0("Outputs/1_DE_infection_modeling/GSEA/",data,"/resultsALL_with_genes.txt"), header=TRUE, sep=" ")
}else {
    results_with_genes <- read.table(paste0("Outputs/1_DE_infection_modeling/GSEA/",data,"/resultsALL_with_genes.txt"), header=TRUE, sep=" ")
    rownames(results_with_genes) <- results_with_genes$peak_ids
}

## connect gene_id with counts matrix as a column.
## add feature id to counts_matrices and results_with_genes from rownames
reads_whole$feature_id <- rownames(reads_whole)
results_with_genes$feature_id <- rownames(results_with_genes)

## merge
res <- full_join(reads_whole, results_with_genes, by="feature_id")

## further filter to only include those genes that are in the immune pathways
res_filt <- res[res$symbol %in% pathway_genes, ]

## remove unecessary columns and only keep reads and symbol
res_filt <- dplyr::select(res_filt, -c(betas, t.statistic, pvalues, fdrs, qvals))
rownames(res_filt) <- res_filt$feature_id

################################################
## read in STR genotypes (used in matrixeqtl) ##
################################################

STRgenotypes = as.data.frame(fread(paste0("Inputs/QTL_mapping/STR_genotypes.txt")))
rownames(STRgenotypes)<-STRgenotypes$V1; STRgenotypes$V1<- NULL
STRgenotypes$str <- rownames(STRgenotypes)

## read in top STRs
## for NI condition
NI_topSTRs <- as.data.frame(fread(paste0("Outputs/3_QTL_mapping/STR-QTL_mapping/",data,"/NI/", data,"_NI/results_best_STRs_with_qval.txt")))
NI_topSTRs$V1<- NULL

## only keep str and gene columns
NI_topSTRs <- dplyr::select(NI_topSTRs, snps, gene)
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
FLU_topSTRs <- as.data.frame(fread(paste0("Outputs/3_QTL_mapping/STR-QTL_mapping/",data,"/Flu/", data,"_Flu/results_best_STRs_with_qval.txt")))
FLU_topSTRs$V1<- NULL

## only keep snp and gene columns
FLU_topSTRs <- dplyr::select(FLU_topSTRs, snps, gene)
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

## now create matrix where you have the top snp for each feature in  each condition, but together in one matrix
both_gene_STRgenotype <-full_join(NI_gene_STRgenotype, FLU_gene_STRgenotype, by="feature")
rownames(both_gene_STRgenotype) <- both_gene_STRgenotype$feature; both_gene_STRgenotype$feature <- NULL

################################################
## read in SNP genotypes (used in matrixeqtl) ##
################################################

SNPgenotypes = as.data.frame(fread(paste0("Inputs/QTL_mapping/SNP_genotypes.txt")))
rownames(SNPgenotypes)<-SNPgenotypes$V1; SNPgenotypes$V1<- NULL
SNPgenotypes$snp <- rownames(SNPgenotypes)

## read in top SNPs for NI condition
NI_topSNPs <- as.data.frame(fread(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/NI/", data,"_NI/results_best_SNPs_with_qval.txt")))
NI_topSNPs$V1<- NULL

## only keep snp and gene columns
NI_topSNPs <- dplyr::select(NI_topSNPs, snps, gene)
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
FLU_topSNPs <- as.data.frame(fread(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/Flu/", data,"_Flu/results_best_SNPs_with_qval.txt")))
FLU_topSNPs$V1<- NULL

## only keep snp and gene columns
FLU_topSNPs <- dplyr::select(FLU_topSNPs, snps, gene)
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

## now create matrix where you have the top snp for each feature in  each condition, but together in one matrix
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
## will remove the gene symbol which we will add back later
res_filt=res_filt[,which(colnames(res_filt) %in% rownames(cols))]
#order the reads count matrix.
res_filt=res_filt[,order(colnames(res_filt))]

#this should be 0. It is a check to ensure that all the metadata and read counts matrix have all the same samples.
length(which(rownames(cols)!=colnames(res_filt)))

## lastly, filter res_filt to make sure that only those features with a SNP and STR are included  (some may not due to window size, etc.)
res_filt <- res_filt[rownames(res_filt) %in% rownames(both_gene_SNPgenotype), ] 
res_filt <- res_filt[rownames(res_filt) %in% rownames(both_gene_STRgenotype), ]

## also filter both_gene_SNPgenotype and both_gene_STRgenotype
both_gene_SNPgenotype_filt <- both_gene_SNPgenotype[rownames(both_gene_SNPgenotype) %in% rownames(res_filt), ] 
both_gene_STRgenotype_filt <- both_gene_STRgenotype[rownames(both_gene_STRgenotype) %in% rownames(res_filt), ] 

##make sure all are same dimension and print number
dim(res_filt)[1]
dim(both_gene_SNPgenotype_filt)[1]
dim(both_gene_STRgenotype_filt)[1]

##################################
#### Add genotype to metadata  ###
##################################

features <- rownames(res_filt)
`%nin%` = Negate(`%in%`)

for(j in 1:length(features)){
  if(j%%100 == 0) print(j)
  feat <- features[j]

  meta_data_loop <- cols
  meta_data_loop$indiv_ID <- rownames(meta_data_loop)

  ## add genotype into meta data
  meta_data_loop$SNP <- as.numeric(both_gene_SNPgenotype_filt[feat, rownames(meta_data_loop)])
  meta_data_loop$STR <- as.numeric(both_gene_STRgenotype_filt[feat, rownames(meta_data_loop)])

  ########################################################################
  ## Get residuals after accounting for top variants for each feature  ###
  ########################################################################

  design=model.matrix(~0 + SNP + STR ,data=meta_data_loop)

  #Remove individuals with an NA genotype
  reads_loop <- res_filt[,colnames(res_filt) %in% rownames(design)]
  fit <-lmFit(reads_loop,design)
  fit <- eBayes(fit)
  ## get residuals
  residuals <- residuals.MArrayLM(object = fit, reads_loop)  

  ## only select row for feature we corrected for
  corrected_matrix_row <- as.data.frame(residuals[rownames(residuals) == feat, ])
  colnames(corrected_matrix_row) <- feat

  if(j==1){
  corrected_matrix <- corrected_matrix_row
  } else {
  corrected_matrix <- as.data.frame(merge(as.data.frame(corrected_matrix), as.data.frame(corrected_matrix_row), by="row.names", all=TRUE))
  rownames(corrected_matrix) <- corrected_matrix$Row.names; corrected_matrix$Row.names <- NULL
  }
}


final_corrected_matrix <- t(corrected_matrix)
dim(final_corrected_matrix)[1]==length(features) 

## fill in missing values with mean for each row/feature
k <- which(is.na(final_corrected_matrix), arr.ind=TRUE)
final_corrected_matrix[k] <- rowMeans(final_corrected_matrix, na.rm=TRUE)[k[,1]]
  
write.table(final_corrected_matrix, paste0(out_dir, "/corrected_matrix_residuals.txt"))
