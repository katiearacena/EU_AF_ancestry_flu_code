## prepare QTL PIP annotation for LDSC
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(doParallel))
options(scipen=999)

args <- commandArgs(trailingOnly = TRUE)
susie_prior <- args[1]

# ## Example settings
# susie_prior <- "dist"

susie_dir <- "/project2/xinhe/kevinluo/flu_regulatoryQTLs/susie"
annot_dir <- "/project2/xinhe/kevinluo/flu_regulatoryQTLs/sldsc_results/annot/annot_bed/"

dir.create(annot_dir, showWarnings = F, recursive = T)

phenotype_list <- c("RNAseq_Flu", "RNAseq_NI",
                    "ATACseq_Flu", "ATACseq_NI",
                    "H3K27ac_Flu", "H3K27ac_NI",
                    "H3K27me3_Flu", "H3K27me3_NI",
                    "H3K4me1_Flu", "H3K4me1_NI",
                    "H3K4me3_Flu", "H3K4me3_NI")

# phenotype_list <-  c(phenotype_list,
#                      paste0("WGBS_Flu_full_E", 3:7), paste0("WGBS_NI_full_E", 3:7))
#
# phenotype_list <- c("WGBS_Flu", "WGBS_NI")

for(phenotype_name in phenotype_list) {

  ## Parameters
  cat("Phenotype name:", phenotype_name, "\n")

  ## QTL annotations
  ### Load Susie fine-mapping results
  susie_output_dir <- paste0(susie_dir, "/output/", phenotype_name, "/susie_", susie_prior, "_prior")
  SusieResult <- readRDS(paste0(susie_output_dir, "/susie_", phenotype_name, "_", susie_prior, "_prior_sumstats.rds"))

  ### take max PIP per SNP
  SusieResult <- SusieResult[order(SusieResult$PIP,decreasing = T),]
  SusieResult <- SusieResult[!duplicated(SusieResult$snps),]
  if(! grepl("chr", SusieResult$chr[1])){
    SusieResult$chr <- paste0("chr", SusieResult$chr)
  }

  SusieResult$pos <- as.numeric(SusieResult$pos)

  ### save continuous PIP as annotations for LDSC
  QTL_continuous_PIP.bed <- data.frame(chr = SusieResult$chr,
                                       start = SusieResult$pos - 1,
                                       end = SusieResult$pos,
                                       PIP = round(SusieResult$PIP,6))
  QTL_continuous_PIP.bed$chr <- factor(QTL_continuous_PIP.bed$chr, levels = paste0("chr", 1:22))
  QTL_continuous_PIP.bed <- QTL_continuous_PIP.bed[order(QTL_continuous_PIP.bed$chr, QTL_continuous_PIP.bed$start), ]
  QTL_continuous_PIP.bed <- unique(QTL_continuous_PIP.bed)
  cat(nrow(QTL_continuous_PIP.bed), "SNPs with PIP annotations. \n")
  fwrite(QTL_continuous_PIP.bed, paste0(annot_dir, "/", phenotype_name, "_PIP_", susie_prior, "_prior_FDR0.1.bed"), sep = "\t")

}

cat("PIP based continuous annotation bed files saved at:", annot_dir, "\n")
