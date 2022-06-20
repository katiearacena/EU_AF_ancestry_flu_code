### To integrate QTL across datatypes and collect ALL snps

## Load required libraries.
library(dplyr)
library(parallel)
library(data.table)

args<-commandArgs(TRUE)

## Declare data type
QTL_type <- as.numeric(args[1])
condition <- args[2]

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
setwd('../../../../')

## Create directory structure to save outputs.
system(paste0("mkdir -p Outputs/",folder,"/SNP-QTL_integration_all_snps/", condition, "/"))
out_dir <- paste0("Outputs/",folder,"/SNP-QTL_integration_all_snps/", condition, "/")

datatypes <- c("RNAseq", "ATACseq", "H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "WGBS")

best_snps<- list()
full_res <- list()

## load best snps from QTL mapping
for (i in 1:length(datatypes)){
  data_i <- datatypes[i]
  print(data_i)
  best_snps[[i]] <- as.data.frame(fread(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/", data_i, "/", condition, "/", data_i, "_", condition, "/results_best_SNPs_with_qval.txt")))
  rownames(best_snps[[i]]) <- best_snps[[i]]$V1; best_snps[[i]]$V1 <- NULL
}

## load full results from QTL mapping
for (i in 1:length(datatypes)){
  data_i <- datatypes[i]
  print(data_i)
  full_res[[i]] <- as.data.frame(fread(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/", data_i, "/", condition, "/", data_i, "_", condition, "/raw_results/result_original.txt")))
  rownames(full_res[[i]]) <- full_res[[i]]$V1; full_res[[i]]$V1 <- NULL

}

### PRINT FOR WHICH DATASET WE ARE STARTING WITH
print(datatypes[QTL_type])


## get eQTL full results
QTL <- full_res[[QTL_type]]
## get qvalue threshold for eQTL. do this by thresholding by qvalue and then selecting the max pvalue that passes that threshold 
sig_eQTL_threshold <- as.vector(summary(best_snps[[QTL_type]][best_snps[[QTL_type]]$qvalues < .10, ]$pvalue)[6])
## make sure you collect those that are less than AND EQUAL to this max pvalue. 
sig_QTL <- QTL[QTL$pvalue <= sig_eQTL_threshold, ]

## make sure the number of genes and the number of eQTL from best snps using qvalue is the same (=TRUE).
length(unique(sig_QTL$gene))== dim(best_snps[[QTL_type]][best_snps[[QTL_type]]$qvalues < .10, ])[1]

## turn genes into a vector
genes <- as.vector(unique(sig_QTL$gene))

print("number of features")
length(genes)
print("number of snps for these features")
dim(sig_QTL)[1]

## create a list which contains dataframes for each gene, that contain all significant snps for that gene
sig_QTL_gene_list <- list()

for (i in 1:length(genes)){
     gene_i <- genes[i]
     snp_genei_pairs <- sig_QTL[sig_QTL$gene == gene_i, ]
     sig_snp_genei_pairs <- snp_genei_pairs[snp_genei_pairs$pvalue <= sig_eQTL_threshold, ]
     sig_snp_genei_pairs_snps <- as.vector(sig_snp_genei_pairs$snps)
     sig_QTL_gene_list[[i]] <- c(sig_snp_genei_pairs_snps)
}
names(sig_QTL_gene_list) <- as.vector(genes)
## now have list with a dataframe for each gene, with all snps that are signifcant eQTL in it.


## only keep rows included in sig_QTL
snps_res <- list()
snps_union <- as.vector(sig_QTL$snps)

for (j in 1:length(full_res)){
    file_name <- full_res[[j]]
    ## only keep rows included in snps_union
    snp_rows <- file_name[file_name$snps %in% snps_union, ]
    ## select best (lowest) pvalue for each snp
    snp_pval_ord <- snp_rows[order(snp_rows[, "snps"], snp_rows[, "pvalue"]),]
    snps_res[[j]] <-snp_pval_ord[!duplicated(snp_pval_ord$snps),]
    }

names(snps_res) <- datatypes
## set function to search for and collect ALLs pvalues for each snp
code1 <- function(my_number) {
	snp_k <- gene_df[my_number]
	snp_rows <- file_name[file_name$snps == snp_k,] 
    
    if(dim(snp_rows)[1]>0){
	## collect top pvalue for snp. 
	    snp_rows_ordered <- snp_rows[order(snp_rows[, "pvalue"]),][1, ]
        top_pval_among_snps_and_snpid <- c(as.character(snp_rows_ordered$snps), snp_rows_ordered$pvalue)
    } else {
        top_pval_among_snps_and_snpid<-c(as.character(snp_k), "NA")
    }
    	return(top_pval_among_snps_and_snpid)
}

## create matrix to store results, columns will be gene pvalues and gene -snps and rows will be datatypes.

## for each gene in sig_QTL_gene_list
for (i in 1:length(sig_QTL_gene_list)){

  gene_df <- sig_QTL_gene_list[[i]]
  if(i%%1==0)print(i)

## create empty vector to store pvalue for the best snp from each datatype in
   data_df <- matrix(ncol=2,nrow=length(snps_res))
## and for each datatype in snps_res
  for (j in 1:length(snps_res)){
    file_name <- snps_res[[j]]

    ## look for each snp in gene_df and find all rows containing it.  
    all_pvals_for_gene <- do.call("rbind",mclapply(1:length(gene_df), FUN = code1))
    ## now have list of top pvalues for each feature with snp.
    all_pvals_for_gene <- as.data.frame(all_pvals_for_gene)
    all_pvals_for_gene$V2 <- as.numeric(levels(all_pvals_for_gene$V2))[all_pvals_for_gene$V2]
    all_pvals_for_gene$V1 <- (as.character(all_pvals_for_gene$V1))
    ## instead of getting top pvalue instead save all pvalues but still order them
    all_pvals_for_gene <- all_pvals_for_gene[order(all_pvals_for_gene[, 2]), ]
    ## rename colnames
    colnames(all_pvals_for_gene)[1] <- "snp"
    colnames(all_pvals_for_gene)[2] <- c(paste0(names(snps_res[j]), "_pval"))
  
    all_pvals_for_gene$feature_id <- names(sig_QTL_gene_list)[i]
  
    if(j==1){
        data_df <- all_pvals_for_gene
    }else{
        data_df <- full_join(data_df, all_pvals_for_gene, by = c("snp", "feature_id")) 
    }
  }
## now save this dataframe for each eGene
  if(i==1){
    gene_data_df <- data_df
  } else {
    gene_data_df <- rbind(gene_data_df, data_df)
  }
}

fwrite(gene_data_df, paste0(out_dir, datatypes[QTL_type], "_", condition, "_QTL.integ.ALL.snps.txt"), na="NA")