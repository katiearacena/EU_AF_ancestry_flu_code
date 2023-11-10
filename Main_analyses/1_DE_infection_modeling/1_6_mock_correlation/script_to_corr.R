#Correlation between mock, NI and Flu. for marks in which they are available. 
#Ideally correlation is strong. Correlation between mock and NI on average is stronger than comparing different individuals.

# Load libraries
library(limma)
library(edgeR)
library(ggplot2)
library(reticulate)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(ggfortify)
library(grid)
library(cowplot)
library(bsseq)
library(stringr)
library(naniar)
library(bsseq)
library(gdata)
library(matrixStats)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(tidyverse)
library(ggridges)
library(PNWColors)


inDir= "/project2/lbarreiro/users/katie/EU_AF_ancestry_flu/"

datatypes <- c("RNAseq","ATACseq","H3K27ac","H3K27me3","H3K4me1", "H3K4me3", "WGBS")
my_colors <- pnw_palette(name="Sunset",n=7,type="discrete")

## batch/age corrected reads do not contain mock since each condition was performed separetly.

for (i in 1:length(datatypes)){
  data_i <- datatypes[i]
  print(data_i)
  
  metadata <- read.table(paste0(inDir, "paper_pipeline/Inputs/metadata/", data_i, "_metadata.txt"))
  
  if(data_i =="RNAseq"){
    
    metadata <- metadata[which(!metadata$Genotyping_ID=="EU122_iPSCs_D53"),]
    metadata <- metadata[which(!metadata$Genotyping_ID=="EU122_iPSCs_D81"),]
  }
  
  if(data_i=="WGBS"){
    load(paste0(inDir, "/WGBS/BSobj.FILTERED_WITH_MOCK.RData"))
    cts <- getMeth(BSobj.fit.nolow2, type="raw")
    reads <- as.data.frame(na.omit(cts))
    colnames(reads) <- sub("AF22", "EU22", colnames(reads))
    colnames(reads) <- sub("AF36", "EU36", colnames(reads))
    colnames(reads) <- sub("AF38", "EU38", colnames(reads))
    reads <- select(reads, -c(EU09_Mock))
  }else if (data_i=="RNAseq"){
    reads <- read.table(paste0(inDir, "paper_pipeline/Inputs/counts_matrices/", data_i, "_filtered.counts.txt"), header=T, row.names=1)
  } else {
    reads <- read.table(paste0(inDir, "mock_counts/", data_i, "_counts_mock.txt"), header=T, row.names=1)
  }

  
  ## Factorize certain columns from meta data
  ## *Batch MUST be a factor so you get multiple estimates!*
  metadata$Batch <- as.factor(metadata$Batch)
  metadata$Genotyping_ID <- as.factor(metadata$Genotyping_ID)
  ## Set NI as first level so that it is used as the reference.
  metadata$Condition = factor(metadata$Condition, levels=c("NI","Flu","Mock"))
  
  ## Set original names in case reordering is necessary.
  rownames(metadata)<- sub("_H3K27ac", "", rownames(metadata))
  rownames(metadata)<- sub("_H3K27me3", "", rownames(metadata))
  rownames(metadata)<- sub("_H3K4me1", "", rownames(metadata))
  rownames(metadata)<- sub("_H3K4me3", "", rownames(metadata))
  metadata$Sample <- sub("_H3K27ac", "", metadata$Sample)
  metadata$Sample<- sub("_H3K27me3", "", metadata$Sample)
  metadata$Sample<- sub("_H3K4me1", "", metadata$Sample)
  metadata$Sample<- sub("_H3K4me3", "", metadata$Sample)
  reorder_names <- rownames(metadata)
  
  ## Read in reads counts and make sure that the dimensions of the reads and metadata match.
  ## Check to see if the length of reorder_names (samples) matches the number of columns in the reads file (samples).
  ## Need number of samples to match. If they match.
  dim(metadata)
  dim(reads)
  length(reorder_names)
  
  ## If these values match, subset reads to include only the samples in metadata file.
  ## This will remove the IPSC samples in my data set. 
  if(length(reorder_names) == dim(reads)[2]){
    reads <- reads[reorder_names]
    ## If the values do not match, remove those rows from the reads counts file (this well remove Mock samples)  
  }else{
    correct_names <- colnames(reads)
    metadata <- metadata[rownames(metadata) %in% correct_names,]
    reorder_names <- rownames(metadata)
    reads <- reads[reorder_names]
  }
  
  #################################
  ## VOOM NORMALIZE READ COUNTS  ##
  #################################
  if(data_i=="WGBS"){
    reads_t <- t(reads)
    
  } else {
    
    dge <- DGEList(counts = reads)
    dge <- calcNormFactors(dge)
    design = model.matrix(~0 + Batch, data = metadata)
    design <- design[, colSums(design != 0) > 0]
    v <- voom(dge, design, plot = FALSE)
    
    reads_t <- t(v$E)
  }
  
  reads_t <- as.data.frame(reads_t)
  
  ##################################
  #### For MOCK V NI comparison ####
  ##################################
  
  ## subset on those for which there is mock data.
  mock_samples <- as.vector(as.character(metadata[metadata$Condition == "Mock", ]$Genotyping_ID))

  if(data_i=="H3K27ac"){
    #metadata_sub <- metadata_sub[rownames(metadata_sub) != "AF14_Mock", ]
    mock_samples <- mock_samples[mock_samples != "AF14"]
   } else if(data_i=="H3K27me3"){
    mock_samples <- mock_samples[mock_samples != "AF14"] 
   }else if(data_i=="H3K4me1"){
    mock_samples <- mock_samples[mock_samples != "AF14"]
   } else if(data_i=="H3K4me3"){
    mock_samples <- mock_samples[mock_samples != "AF14"]
   }

  ## subset all those with these sample id
  metadata_mock <- metadata[metadata$Genotyping_ID %in% mock_samples, ]
  
  ## subset reads to obtain only these for NI and mock
  reads_mock <- reads_t[rownames(reads_t) %in% metadata_mock$Sample, ]
  
  ## only keep those with NI or Mock
  conditions <- c("NI", "Mock")
  reads_mock <- reads_mock[grep(paste(conditions, collapse="|"), rownames(reads_mock)), ]
  
  reads_mock_t <- as.data.frame(t(reads_mock))
  
  ## correlate
  corr_Mock_plots_temp <- list()
  
  for (j in 1:length(mock_samples)){

    ind_i <-mock_samples[j] 
    print(ind_i)
      
    ## subset this individual from the cpg matrix .
    
    ind_cor <- reads_mock_t[grepl(paste0(ind_i), colnames(reads_mock_t))]
    ind_cor[, 1] <- as.numeric(as.character(ind_cor[, 1]))
    ind_cor[, 2] <- as.numeric(as.character(ind_cor[, 2]))
    
    cor <- cor.test(ind_cor[, paste0(ind_i, "_Mock")], ind_cor[, paste0(ind_i, "_NI")], method = c("pearson"))
    cor_info<- cbind(ind_i, cor$estimate, cor$p.value)
      

    if(j==1){
      cor_mock_df <- cor_info
    } else {
      cor_mock_df <- rbind(cor_mock_df, cor_info)
    }
  }
  
  cor_mock_df <- as.data.frame(cor_mock_df)
  colnames(cor_mock_df) <- c("sample", "pearsons_corr_coef", "pvalue")
  cor_mock_df$datatype <- paste0(data_i)
  cor_mock_df$corr <- "Mock_v_NI"
  
  ##################################
  #### For NI V Flu comparison  ####
  ##################################
  
  Flu_samples <- as.vector(as.character(metadata[metadata$Condition == "Flu", ]$Genotyping_ID))
  
  ## only keep those with NI or Mock
  conditions <- c("NI", "Flu")
  reads_NI_Flu <- reads_t[grep(paste(conditions, collapse="|"), rownames(reads_t)), ]
  reads_NI_Flu_t <- as.data.frame(t(reads_NI_Flu))
  colnames(reads_NI_Flu_t) <- gsub("_dup", "", colnames(reads_NI_Flu_t))
  
  ## correlate
  corr_Flu_plots_temp <- list()
  
  for (j in 1:length(Flu_samples)){
    
    ind_i <-Flu_samples[j] 
    print(ind_i)
    
    ## subset this individual from the cpg matrix .
    
    ind_cor <- reads_NI_Flu_t[grepl(paste0(ind_i), colnames(reads_NI_Flu_t))]
    ind_cor[, 1] <- as.numeric(as.character(ind_cor[, 1]))
    ind_cor[, 2] <- as.numeric(as.character(ind_cor[, 2]))
    
    cor <- cor.test(ind_cor[, paste0(ind_i, "_Flu")], ind_cor[, paste0(ind_i, "_NI")], method = c("pearson"))
    cor_info<- cbind(ind_i, cor$estimate, cor$p.value)
    
    if(j==1){
      cor_Flu_df <- cor_info
    } else {
      cor_Flu_df <- rbind(cor_Flu_df, cor_info)
    }
  }
    
  cor_Flu_df <- as.data.frame(cor_Flu_df)
  colnames(cor_Flu_df) <- c("sample", "pearsons_corr_coef", "pvalue")
  cor_Flu_df$datatype <- paste0(data_i)
  cor_Flu_df$corr <- "Mock_v_Flu"
  
  if(i==1){
    cor_df_all <- rbind(cor_mock_df, cor_Flu_df)
  } else {
    cor_df_all <- rbind(cor_df_all, cor_mock_df, cor_Flu_df)
  }  

}

## make sure all the correct classes
cor_df_all$pearsons_corr_coef <- as.numeric(as.character(cor_df_all$pearsons_corr_coef))
cor_df_all$pvalue <- as.numeric(as.character(cor_df_all$pvalue))
cor_df_all$datatype <- factor(cor_df_all$datatype, levels=c("RNAseq","ATACseq","H3K27ac","H3K27me3","H3K4me1", "H3K4me3", "WGBS"))
cor_df_all$corr <- factor(cor_df_all$corr, levels=c("Mock_v_NI", "Mock_v_Flu"))


write.table(cor_df_all, file=paste0("/project2/lbarreiro/users/katie/EU_AF_ancestry_flu/mock_counts/corr/all_data_corr_coef_df.txt"))