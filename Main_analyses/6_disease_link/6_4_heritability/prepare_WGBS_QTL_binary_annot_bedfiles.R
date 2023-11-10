## prepare regulatory QTL annotations in BED format for S-LDSC analysis
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
thresh_FDR <- as.numeric(args[1])
outdir <- args[2]

QTL_res_dir <- "/project2/xinhe/kevinluo/flu_regulatoryQTLs/QTL_results/"

cat("thresh_FDR =", thresh_FDR, "\n")

annot_list <- paste0(rep(paste0("WGBS_", c("NI","Flu")),4), rep(paste0("_full_E", 3:7), each = 2))

for(annot_name in annot_list) {
  cat("loading", annot_name, '\n')
  # Use SNP level FDR < 10% cutoff
  sig_annot <- fread(paste0(QTL_res_dir, "/subset_WGBS/", annot_name, "_sig_FDR", thresh_FDR, ".txt.gz"))

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
