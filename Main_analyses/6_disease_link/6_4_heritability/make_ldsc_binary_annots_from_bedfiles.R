## make binary annotation file with multiple annotation bed files and bim file from Plink
library(optparse)
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))

# Process the command-line arguments --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
parser <- OptionParser()
parser <- add_option(parser,"--beddir",type="character",default=NULL)
parser <- add_option(parser,"--bimfile",type="character",default=NULL)
parser <- add_option(parser,"--out",type="character", default="ldsc.annot.gz")
out    <- parse_args(parser)

beddir <- out$beddir
bimfile <- out$bimfile
outfile <- out$out

outdir <- dirname(outfile)
if(!dir.exists(outdir)) dir.create(outdir, showWarnings = F, recursive = T)

cat("Make LDSC annotations from bed files...\n")
cat("bim file:", bimfile, "\n")

bedfiles <- list.files(path = beddir, pattern = "\\.bed$", full.names = TRUE)
cat(length(bedfiles), "bed files. \n")

bim_data <- fread(bimfile, col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
annot_ldsc <- bim_data[, c("CHR", "BP", "SNP", "CM")]
annot_ldsc.gr <- makeGRangesFromDataFrame(annot_ldsc, seqnames.field = "CHR", start.field = "BP", end.field = "BP", keep.extra.columns = T)

annot_ldsc$base <- 1

for(bedfile in bedfiles){
  annot_name <- gsub(".bed", "", basename(bedfile))
  cat("Make annot for:", annot_name, "\n")
  annot_bed <- fread(bedfile, select = c(1:3), col.names = c("chr", "start", "end"))
  annot_bed$chr <- gsub("chr", "", annot_bed$chr)
  annot_bed$start <- as.integer(annot_bed$start) + 1
  annot_bed$end <- as.integer(annot_bed$end)
  annot_bed.gr <- makeGRangesFromDataFrame(annot_bed, keep.extra.columns = T)
  overlaps <- findOverlaps(query = annot_ldsc.gr, subject = annot_bed.gr)
  annot_ldsc[, annot_name] <- 0
  annot_ldsc[queryHits(overlaps), annot_name] <- 1
}

fwrite(annot_ldsc, outfile, sep = "\t")
cat("Done. LDSC annotation file saved at: ", outfile, "\n")

