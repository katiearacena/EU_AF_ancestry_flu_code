library(stringi)
library(stringr)
library(data.table)
library(dplyr)
set.seed(2021)

args<-commandArgs(TRUE)

## Declare data type
QTLnum <- args[1]
condition <- args[2]

if(length(args)==0)
  {
    print("WARNING: No arguments supplied.")
  }
print(args)

print(QTLnum)
print(condition)

#############################################################
## CREATE DIRECTORY STRUCTURE, SET FUNCTIONS, & LOAD DATA  ##
#############################################################

folder = "3_QTL_mapping"
## Set directory 3 steps above script.
setwd('../../../')

## Create directory structure to save outputs.
system(paste0("mkdir -p Outputs/",folder,"/QTL_integ_TF_enrichments/"))
out_dir <- paste0("Outputs/",folder,"/QTL_integ_TF_enrichments/")

## load TF files
footprints <- as.data.frame(fread(paste0("Inputs/footprints/",condition, "_footprints_motifs_overlapped.txt")))
colnames(footprints) <- c("chr", "start", "end", "TF_motif")

chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")

## load QTL integration union lists
QTL_union_vec <- read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_integration_all_snps/", condition, "/union_lists_for_TF/", QTLnum, "_", condition, "_union.txt"))
colnames(QTL_union_vec)  <- c("union_snps")
## remove feature ids from the snps
snps_from_union_vec <- as.data.frame(stringr::str_extract(as.vector(QTL_union_vec$union_snps), "[^_]*_[^_]*_[^_]*_[^_]*"))
colnames(snps_from_union_vec) <- c("snps")
print(paste0("number of snps for ",QTLnum, " = ", dim(snps_from_union_vec)[1]))
##only keep unique snps because same snp could be QTL for different features
uniq_snps_from_union_vec <- as.data.frame(unique(snps_from_union_vec$snp))
colnames(uniq_snps_from_union_vec) <- c("snps")
print(paste0("number of unique snps for ",QTLnum, " = ", dim(uniq_snps_from_union_vec)[1]))

## obtain snp location for QTL results
uniq_snps_from_union_vec$chr <- str_extract(uniq_snps_from_union_vec$snps, "[^_]+")
uniq_snps_from_union_vec$chr <- paste0("chr", uniq_snps_from_union_vec$chr)
uniq_snps_from_union_vec$pos <- as.numeric(sub("^[^_]*_([^_]*).*", "\\1", uniq_snps_from_union_vec$snps))

## set function to collect if snp falls within a TF footprint or not.
get_region <- function(vec, id) {
if(length(.footprints <- which(vec >= footprints$start & vec <= footprints$end & id == footprints$chr)[1])) .footprints else NA
}

sig_QTL <- uniq_snps_from_union_vec
    
## set empty dataframe to store results
sig_QTL_with_overlap <- list()

for (chr in 1:length(chromosomes)){
    chr_i<- chromosomes[chr]
    print(chr_i)
    chr_i_QTL <- sig_QTL[sig_QTL$chr == chr_i, ]

    if(dim(chr_i_QTL)[1] > 0){
        chr_i_QTL$overlapping_footprint <- as.character(footprints$TF_motif[mapply(get_region, chr_i_QTL$pos, chr_i)])
    }
    sig_QTL_with_overlap[[chr]] <- chr_i_QTL
}
sig_QTL_with_overlap <- do.call(rbind, sig_QTL_with_overlap)

sig_QTL_with_overlap$overlapping_footprint[is.na(sig_QTL_with_overlap$overlapping_footprint)] <- "no_overlap"

print(dim(sig_QTL_with_overlap)[1]==dim(sig_QTL)[1])

## write file
saveRDS(sig_QTL_with_overlap,  file = paste0(out_dir, condition, "/", QTLnum, ".TF.overlaps.RDS"))