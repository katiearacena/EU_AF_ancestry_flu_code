
library(tidyverse)
library(ggplot2)
library(dplyr)



## load from cluster (midway3)
input_dir <- "/Volumes/project/lbarreiro/USERS/katie/EU_AF_ancestry_flu/local_ancestry/rfmix/outputs/"
out_dir <- "/Volumes/project/lbarreiro/USERS/katie/EU_AF_ancestry_flu/local_ancestry/rfmix/LA_calls_for_popDE/"
## load but remove first line which details subpopulation order and codes
## CEU=0
## YRI=1

sample_file <- read.delim(paste0(input_dir, "ref_1.msp.tsv"), header=TRUE, skip=1)
samples_d1 <- colnames(sample_file[7:ncol(sample_file)])
samples_d2 <- str_subset(samples_d1, "1.0")
samples_list <- gsub(".0", "", samples_d2, fixed=TRUE)


chrs <- c(1:22)
LA_calls_samples_list <- list()
sample_names_vec <- vector()


for (s in 1:length(samples_list)){
  sample_i <- samples_list[s]
  
  for (i in 1:length(chrs)){
    
    chr_i <- chrs[i]
    queery_samples <- read.delim(paste0(input_dir, "ref_", chr_i, ".msp.tsv"), header=TRUE, skip=1)
    gnomix <- queery_samples %>% select(X.chm, spos, epos, starts_with(as.character(sample_i)))
    
    sample_i2 <- gsub('Epi_', "", sample_i)
    sample_cut <- stringr::str_extract(sample_i2, "[^_]*_[^_]*")
    sample_names_vec[s] <- sample_cut
    
    ## get sum of rows. will be either 0= 2 copies CEU, 1, or 2= 2 copies YRI ancestry at locus
    gnomix2 <- gnomix %>% mutate(LA_call = select(., starts_with(as.character(sample_i))) %>% rowSums(na.rm = TRUE))
    
    ## make 1 file with all chrs for each individual
    if(i==1){
      chr_df <- gnomix2
    } else {
      chr_df <- rbind(chr_df, gnomix2)
    }
  }
  
  chr_df_filt <- select(chr_df, "X.chm", "spos", "epos", "LA_call")
  LA_calls_samples_list[[s]] <- chr_df_filt
  
  #write_tsv(chr_df, file = paste0(out_dir, sample_cut, "_LA_calls_for_popDE_ALL_chrs_ref_rfmix_tagore.bed"), col_names = F)
}

names(LA_calls_samples_list) <- sample_names_vec 

## first perform some QC checks. 
## are all regions the same across samples? they should be but just check
## see if all other dfs match first dfs

LA_calls_samples_REF <- LA_calls_samples_list[[1]]

for (n in 2:length(LA_calls_samples_list)){
  print(n)
  
  LA_calls_samples_TO_COMPARE = LA_calls_samples_list[[n]]
  
  print(all(LA_calls_samples_REF$spos == LA_calls_samples_TO_COMPARE$spos))
  print(all(LA_calls_samples_REF$epos == LA_calls_samples_TO_COMPARE$epos))
  ## should all be TRUE. an issue if false
}

## now extract relevant columns and make a matrix.
for (n in 1:length(LA_calls_samples_list)){

  LA_calls_sample = LA_calls_samples_list[[n]]
  colnames(LA_calls_sample)[4] = names(LA_calls_samples_list)[n]
  if(n==1){
    df <- LA_calls_sample
  } else {
    LA_call_only <- select(LA_calls_sample, names(LA_calls_samples_list)[n])
    df <- cbind(df, LA_call_only)
  }
}

colnames(df)[1] <- "chr"
##save
write.table(df, file=paste0(out_dir, "LA_matrix_for_popDE.txt"))

