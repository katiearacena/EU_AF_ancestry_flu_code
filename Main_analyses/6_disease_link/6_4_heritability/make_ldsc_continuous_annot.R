## make continuous annotation file from bim file and annotation bed file

# Load packages
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

# Process the command-line arguments.
option_list <- list(
  make_option("--bed", action="store", default=NA, type='character',
              help="bed file of the annotation"),
  make_option("--bim", action="store", default=NA, type='character',
              help="1000G bim file "),
  make_option("--out", action="store", default=NA, type='character',
              help="Output file"),
  make_option("--out_format", action="store", default="full-annot", type='character',
              help="Output annotation format. Option: full-annot or thin-annot")
)

opt <- parse_args(OptionParser(option_list=option_list))
bed_file      <- opt$bed
bim_file      <- opt$bim
out_file      <- opt$out
out_format    <- opt$out_format

annot_bed <- fread(bed_file, select = c(1:4), col.names = c("chr", "start", "end", "value"))
annot_bed$chr <- gsub("chr", "", annot_bed$chr)
annot_bed$start <- as.numeric(annot_bed$start) + 1
annot_bed$end <- as.numeric(annot_bed$end)

cat("Make LDSC continuous annot for: ", bed_file, "\n")

bim_data <- fread(bim_file, col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))

annot_ldsc <- bim_data[, c("CHR", "BP", "SNP", "CM")]

annot_ldsc.gr <- makeGRangesFromDataFrame(annot_ldsc, seqnames.field = "CHR", start.field = "BP", end.field = "BP", keep.extra.columns = T)

annot_bed.gr <- makeGRangesFromDataFrame(annot_bed, keep.extra.columns = T)

idx_overlaps <- data.frame(findOverlaps(query = annot_ldsc.gr, subject = annot_bed.gr))

annot_ldsc[, "ANNOT"] <- 0
annot_ldsc[idx_overlaps$queryHits, "ANNOT"] <- annot_bed$value[idx_overlaps$subjectHits]

dir.create(dirname(out_file), showWarnings = F, recursive = T)

if(grepl("thin", out_format)){
  fwrite(data.frame(ANNOT = annot_ldsc[, "ANNOT"]), out_file, sep = "\t")
}else if(grepl("full", out_format)){
  fwrite(annot_ldsc, out_file, sep = "\t")
}

cat("LDSC continuous annotation file saved at: ", out_file, "\n")

