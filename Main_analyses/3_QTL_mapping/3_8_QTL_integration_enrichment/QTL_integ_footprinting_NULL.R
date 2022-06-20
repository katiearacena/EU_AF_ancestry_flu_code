library(stringi)
library(stringr)
library(data.table)
library(dplyr)
set.seed(2021)

args<-commandArgs(TRUE)

chrnum <- as.numeric(args[1])
condition <- args[2]

if(length(args)==0)
  {
    print("WARNING: No arguments supplied.")
  }
print(args)

print(chrnum)
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

## set function to collect if snp falls within a TF footprint or not.
get_region <- function(vec, id) {
if(length(.footprints <- which(vec >= footprints$start & vec <= footprints$end & id == footprints$chr)[1])) .footprints else NA
}

#####################################################
## Null = any snps that are tested for QTL mapping ## 
#####################################################

## load file which has union of all snps tested for mapping
## load full results from QTL mapping

full_res <- list()

datatypes <- c("RNAseq", "ATACseq", "H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "WGBS")

for (i in 1:length(datatypes)){
  data_i <- datatypes[i]
  print(data_i)
  full_res[[i]] <- as.data.frame(fread(paste0("/Outputs/3_QTL_mapping/SNP-QTL_mapping/", data_i, "/", condition, "/", data_i, "_", condition, "/raw_results/result_original.txt")))
  rownames(full_res[[i]]) <- full_res[[i]]$V1; full_res[[i]]$V1 <- NULL
}

## only select snp ids
all_matrixeqtl_snps <- lapply(full_res, function(x) { dplyr::select(x, snps)})
all_matrixeqtl_snps_df <- do.call(rbind, all_matrixeqtl_snps)
dim(all_matrixeqtl_snps_df)[1]
unique_all_matrixeqtl_snps_df <- unique(all_matrixeqtl_snps_df$snps)
length(unique_all_matrixeqtl_snps_df)

saveRDS(unique_all_matrixeqtl_snps_df, paste0(out_dir, "null_snps_vector.RDS"))

null2_snps_vector <- readRDS(paste0(out_dir, "null_snps_vector.RDS"))
print(paste0("number of unique snps for null#2 for ",condition, " = ", length(null2_snps_vector)))

null2_df <- as.data.frame(null2_snps_vector)
colnames(null2_df) <- c("snps")

## obtain snp location
null2_df$chr <- str_extract(null2_df$snps, "[^_]+")
null2_df$chr <- paste0("chr", null2_df$chr)
null2_df$pos <- as.numeric(sub("^[^_]*_([^_]*).*", "\\1", null2_df$snps))

chr_i<- chromosomes[chrnum]
print(chr_i)
chr_i_QTL <- null2_df[null2_df$chr == chr_i, ]

print(paste0("number of unique snps for null#2 for chr",chrnum, " = ", dim(chr_i_QTL)[1]))

if(dim(chr_i_QTL)[1] > 0){
  chr_i_QTL$overlapping_footprint <- as.character(footprints$TF_motif[mapply(get_region, chr_i_QTL$pos, chr_i)])
}

chr_i_QTL$overlapping_footprint[is.na(chr_i_QTL$overlapping_footprint)] <- "no_overlap"

## write file
saveRDS(chr_i_QTL,  file = paste0(out_dir, condition, "/NULL.chr", chrnum, ".RDS"))