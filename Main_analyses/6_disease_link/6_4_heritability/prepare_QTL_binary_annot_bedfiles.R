## prepare regulatory QTL annotations in BED format for S-LDSC analysis
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
thresh_FDR <- as.numeric(args[1])
outdir <- args[2]

QTL_res_dir <- "/project2/xinhe/kevinluo/flu_regulatoryQTLs/QTL_results/"

cat("thresh_FDR =", thresh_FDR, "\n")

annot_list <- c("ATACseq_Flu", "ATACseq_NI",
                "H3K27ac_Flu", "H3K27ac_NI",
                "H3K27me3_Flu", "H3K27me3_NI",
                "H3K4me1_Flu", "H3K4me1_NI",
                "H3K4me3_Flu", "H3K4me3_NI",
                "RNAseq_Flu", "RNAseq_NI",
                "WGBS_Flu", "WGBS_NI")

for(annot_name in annot_list) {
  cat("loading", annot_name, '\n')
  annot_file <- paste0(QTL_res_dir, "/", annot_name, "/result_original_rsIDs.txt")
  annot_data <- fread(annot_file, drop = 1)

  # Use SNP level FDR < 10% cutoff
  sig_annot  <- annot_data[annot_data$FDR < thresh_FDR, ]
  fwrite(sig_annot, paste0(QTL_res_dir, "/", annot_name, "/sig_result_FDR", thresh_FDR, ".txt.gz"), sep = "\t")
  sig_annot <- fread(paste0(QTL_res_dir, "/", annot_name, "/sig_result_FDR", thresh_FDR, ".txt.gz"), sep = "\t")

  sig_annot <- sig_annot %>%
    separate(snps, c("chr", "pos", "Allele1", "Allele2"), sep = "_", remove = FALSE) %>%
    mutate(chr = paste0("chr", chr))

  sig_annot.bed <- data.frame(chr = sig_annot$chr, start = as.numeric(sig_annot$pos) - 1, end = as.numeric(sig_annot$pos), name = sig_annot$snps)
  sig_annot.bed$chr <- factor(sig_annot.bed$chr, levels = paste0("chr", 1:22))
  sig_annot.bed <- sig_annot.bed[order(sig_annot.bed$chr, sig_annot.bed$start), ]
  sig_annot.bed <- unique(sig_annot.bed)
  cat(nrow(sig_annot.bed), "SNPs in", annot_name, "binary annotation. \n")

  dir.create(outdir, showWarnings = F, recursive = T)
  fwrite(sig_annot.bed, paste0(outdir, "/", annot_name, "_QTLs_FDR", thresh_FDR, ".bed"), sep = "\t")
}


cat("annotation bed files saved at:", outdir, "\n")
