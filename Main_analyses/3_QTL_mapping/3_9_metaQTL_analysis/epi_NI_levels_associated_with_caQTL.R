library(data.table)
library(dplyr)
library(matrixStats)
library(stringr)

args<-commandArgs(TRUE)
data_i <- args[1]

if(length(args)==0)
{
  print("WARNING: No arguments supplied.")
}
print(args)

#############################################################
## CREATE DIRECTORY STRUCTURE, SET FUNCTIONS, & LOAD DATA  ##
#############################################################

folder = "3_QTL_mapping"
## Set directory 3 steps above script.
setwd('../../../')

## Create directory structure to save outputs.
system(paste0("mkdir -p Outputs/",folder,"/metaQTL/"))
out_dir <- paste0("Outputs/",folder,"/metaQTL/")

######################
### Load all files ###
######################

## load epigenetic NI matrixeqtl results
epi_NI_best <- as.data.frame(fread(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/ATACseq/NI/ATACseq_NI/results_best_SNPs_with_qval.txt")))
epi_NI_best$V1 <- NULL

## load SNP genotypes
SNPgenotypes <- as.data.frame(fread(paste0("Inputs/QTL_mapping/SNP_genotypes.txt")))
rownames(SNPgenotypes) <- SNPgenotypes$V1; SNPgenotypes$V1 <- NULL

## epigenetic NI reads
epi_reads <- as.data.frame(fread(paste0("/Outputs/2_ancestry_effects_modeling/batch_corrected_cts_matrices/", data_i, "_batch.age.corrected.txt")))
rownames(epi_reads) <- epi_reads$V1; epi_reads$V1 <- NULL
epi_reads_NI <- epi_reads[ , grepl( "NI" , colnames(epi_reads))]
colnames(epi_reads_NI) <- str_remove(colnames(epi_reads_NI), "_NI")

## load epigenetic feature QTL NI results (full results)
epi_QTL_NI_results <- as.data.frame(fread(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data_i, "/NI/", data_i, "_NI/raw_results/result_original.txt")))
epi_QTL_NI_results$V1 <- NULL

##############################
### create qqnorm function ###
##############################

qqnorm_data <- function(cts) {
quantile_expression <- matrix(, nrow = nrow(cts), ncol = ncol(cts))
for (j in 1:nrow(cts)){
    exp <- cts[j,]
    exp_QN <- qqnorm(exp, plot.it =F)$x
    quantile_expression[j,] <- exp_QN
}
colnames(quantile_expression) <- colnames(cts)
rownames(quantile_expression) <- rownames(cts)
return(quantile_expression)
}

#############################
### apply qqnorm function ###
#############################

epi_qqnorm_FC <- qqnorm_data(epi_reads_NI)

## make sure mean ~0 and variance ~1
rowMeans(epi_qqnorm_FC[1:10, ])
rowVars(epi_qqnorm_FC[1:10, ])

############################################
### make function to create metagenotype ###
############################################

get_meta_geno_QTL<- function(row, qqnorm_r, genotypes){

## first match gene and peak id
colnames(row)[2] <- "feature"

## get snp-epi features that have same snp as snp as test
epi_QTL_NI_results_i <- epi_QTL_NI_results[epi_QTL_NI_results$snps %in% row$snps, ]
## sort based on pvalues
epi_QTL_NI_results_i <- epi_QTL_NI_results_i[order(epi_QTL_NI_results_i$pvalue), ]
## collect feature that has the lowest pvalue
epi_feature <- epi_QTL_NI_results_i[1, "gene"]

row_info <- cbind(row, epi_feature)

## create meta genotype (homo low, het, homo high)

geno_NI <- genotypes[rownames(genotypes) %in% row_info$snp, ]

## convert -9 to NA
geno_NI[geno_NI == -9] <- "NA"
if (row_info$beta > 0){
    geno_NI[geno_NI == 0] <- "l/l"
    geno_NI[geno_NI == 1] <- "l/h"
    geno_NI[geno_NI == 2] <- "h/h"
} else if (row_info$beta < 0){
    geno_NI[geno_NI == 0] <- "h/h"
    geno_NI[geno_NI == 1] <- "l/h"
    geno_NI[geno_NI == 2] <- "l/l"
}

geno_NI_t <- as.data.frame(t(geno_NI))
colnames(geno_NI_t) <- c("geno")
  ## if epi_feature is NA that means that there is not  data for the closest peak so skip
  if(is.na(row_info$epi_feature)==TRUE){
    geno_NI_t$value <- "NA"
    df_to_plot <- geno_NI_t
    print("no peak info")
    
    df_to_plot$epi_feature <- row_info$epi_feature
    rownames(df_to_plot) <- paste0(rownames(df_to_plot), "_NI")
    
  } else {

    for(m in 1:nrow(row_info)){
    NI_peak_row <- qqnorm_r[rownames(qqnorm_r) == row_info[m, ]$epi_feature, ]
    NI_reads <- reshape2::melt(NI_peak_row)
    
    ## merge geno with qqnorm epi read out
    geno_NI_t$variable <- paste0(rownames(geno_NI_t), "_NI")
    NI_reads$variable <- paste0(rownames(NI_reads), "_NI")
    df_to_plot <- merge(geno_NI_t, NI_reads, by="variable")
    
    df_to_plot$epi_feature <- row_info[m, ]$epi_feature
    
    if(m==1){
        df_bind <- df_to_plot
    } else {
        df_bind <- rbind(df_bind, df_to_plot)
    }
    df_to_plot <- df_bind
  }
 }
  df_to_plot$geno <- factor(df_to_plot$geno)

  ## add peak info
  df_to_plot$snp <- row_info$snps
  df_to_plot$ATACseq_peak <- row_info$feature
  df_to_plot$gene_symbol <- row_info$gene_symbol

  return(df_to_plot)
}

###############################
### Subset on q<.10 epi QTL ###
###############################

epi_NI_best_sig <- epi_NI_best[epi_NI_best$qvalues <.10, ]
dim(epi_NI_best_sig)[1]

##############################################
### Apply metagenotype function to epi QTL ###
##############################################

binned_geno_list <- list()
for (i in 1:dim(epi_NI_best_sig)[1]){
    print(i)
    snp_i <- epi_NI_best_sig[i, ]
    binned_geno_list[[i]] <- get_meta_geno_QTL(snp_i, epi_qqnorm_FC, SNPgenotypes)
}

###################
### Save output ###
###################

saveRDS(binned_geno_list, paste0(out_dir, data_i, "_NI_at_caQTL.RDS"))