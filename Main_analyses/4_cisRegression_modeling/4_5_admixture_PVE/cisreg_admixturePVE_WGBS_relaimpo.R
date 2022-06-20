

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(qvalue)
library(relaimpo)

data = "WGBS"

#################################
## Set directory and load data ##
#################################

folder = "4_cisregression_modeling"
## Set directory 3 steps above script.
setwd('../../../')


## Set directory.
out_dir <- paste0("Outputs/",folder,"/cisregSNPsSTRs/", data,"/")

load(file=paste0(out_dir,"/SNP_STR_admixturePVE_WGBS.RData"))


##################################
#### Add genotype to metadata  ###
##################################

both_variants_genes <- intersect(rownames(FLU_STRgenotypes_subset), rownames(FLU_SNPgenotypes_subset))
genes <- both_variants_genes[which(both_variants_genes %in% sig_in_either)]
#genes <- both_variants_genes

#### Determine subset to be analyzed ####
#total_num <- length(genes)
#sub_size <- ceiling(total_num/16)

#start <- (group-1)*sub_size + 1
#end <- group*sub_size

#if(end>total_num){end <- total_num}


#results_FLU = data.frame(matrix(nrow = end-start+1, ncol = 4))
#rownames(results_FLU) <- genes[start:end]
results_FLU = data.frame(matrix(nrow = length(genes), ncol = 4))
rownames(results_FLU) <- genes
colnames(results_FLU) <- c("no_variant", "SNP", "STR", "SNP_STR")

#results_NI = data.frame(matrix(nrow = end-start+1, ncol = 4))
#rownames(results_NI) <- genes[start:end]
results_NI = data.frame(matrix(nrow = length(genes), ncol = 4))
rownames(results_NI) <- genes
colnames(results_NI) <- c("no_variant", "SNP", "STR", "SNP_STR")


for(j in 1:length(genes)){
  if(j%%10 == 0) print(j)
  gene <- genes[j]
  
  meta_data_loop <- cols
  rownames(meta_data_loop) <- sub("_dup", "", rownames(meta_data_loop))
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
