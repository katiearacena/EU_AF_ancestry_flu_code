library(stringi)
library(stringr)
library(data.table)
library(dplyr)
set.seed(2021)

args<-commandArgs(TRUE)

## Declare data type
data_type_i <- args[1]
condition <- args[2]

if(length(args)==0)
  {
    print("WARNING: No arguments supplied.")
  }
print(args)

print(data_type_i)
print(condition)


#############################################################
## CREATE DIRECTORY STRUCTURE, SET FUNCTIONS, & LOAD DATA  ##
#############################################################

folder = "3_QTL_mapping"
## Set directory 3 steps above script.
setwd('../../../')

## Create directory structure to save outputs.
system(paste0("mkdir -p Outputs/",folder,"/condition_specific_QTL_TF_enrichments/"))
out_dir <- paste0("Outputs/",folder,"/condition_specific_QTL_TF_enrichments/")

## load TF files
footprints <- as.data.frame(fread(paste0("Inputs/footprints/",condition, "_footprints_motifs_overlapped.txt")))

colnames(footprints) <- c("chr", "start", "end", "TF_motif")


## load cluster info
cluster <- read.table("Inputs/footprints/clusters.txt")
colnames(cluster) <- c("cluster")
cluster$TF_motif <- rownames(cluster)


## merge TF footprinting with cluster info
footprints_cluster <- full_join(footprints, cluster, by="TF_motif")
footprints_cluster$chr <- as.character(footprints_cluster$chr)

chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")

## load matrix eqtl results. 
NI_q <- as.data.frame(fread(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/", data_type_i,"/NI/", data_type_i, "_NI/results_best_SNPs_with_qval.txt")))
rownames(NI_q) <- NI_q$V1; NI_q$V1 <- NULL
FLU_q <- as.data.frame(fread(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/", data_type_i,"/Flu/", data_type_i, "_Flu/results_best_SNPs_with_qval.txt")))
rownames(FLU_q) <- FLU_q$V1; FLU_q$V1 <- NULL


if(data_type_i =="WGBS"){

    pt1 <- read.table("Outputs/3_QTL_mapping/SNP-QTL_mash/WGBS/part1_outputs/lfsr_output.txt")
    pt2 <- read.table("Outputs/3_QTL_mapping/SNP-QTL_mash/WGBS/part2_outputs/lfsr_output.txt")
    pt3 <- read.table("Outputs/3_QTL_mapping/SNP-QTL_mash/WGBS/part3_outputs/lfsr_output.txt")
    pt4 <- read.table("Outputs/3_QTL_mapping/SNP-QTL_mash/WGBS/part4_outputs/lfsr_output.txt")
    pt5 <- read.table("Outputs/3_QTL_mapping/SNP-QTL_mash/WGBS/part5_outputs/lfsr_output.txt")
    pt6 <- read.table("Outputs/3_QTL_mapping/SNP-QTL_mash/WGBS/part6_outputs/lfsr_output.txt")
    pt7 <- read.table("Outputs/3_QTL_mapping/SNP-QTL_mash/WGBS/part7_outputs/lfsr_output.txt")
    pt8 <- read.table("Outputs/3_QTL_mapping/SNP-QTL_mash/WGBS/part8_outputs/lfsr_output.txt")

    mash <- rbind(pt1, pt2, pt3, pt4, pt5, pt6)
    colnames(mash) <- c("lfsr_FLU", "lfsr_NI")

} else {
    mash <- as.data.frame(fread(paste0("Outputs/3_QTL_mapping/SNP-QTL_mash/", data_type_i, "/results_BOTH.mash.txt")))
    rownames(mash) <- mash$V1; mash$V1 <- NULL
}


NI_q$feature <- NI_q$gene
FLU_q$feature <- FLU_q$gene
mash$feature <- rownames(mash)

NI_mash <- select(mash, feature, lfsr_NI)
FLU_mash <- select(mash, feature, lfsr_FLU)

NI_results <- full_join(NI_q, NI_mash, by="feature")
rownames(NI_results) <- NI_results$feature
FLU_results <- full_join(FLU_q, FLU_mash, by="feature")
rownames(FLU_results) <- FLU_results$feature

### for NI specific
sig_NI_results <- dplyr::filter(NI_results, NI_results$qval < .10)
##NI-infected specific: fdr<0.1 in flu infected & fdr>0.1 in NI condition
R_NI2 <- (dplyr::filter(sig_NI_results, sig_NI_results$lfsr < .10))
R_Flu2 <- (dplyr::filter(FLU_results, FLU_results$lfsr > .10 ))
NI_specific <- intersect(R_Flu2$feature, R_NI2$feature)
#Number of NI specific.
print(length(NI_specific))

### for Flu specific
sig_FLU_results <- dplyr::filter(FLU_results, FLU_results$qval < .10)
##Flu-infected specific: fdr<0.1 in flu infected & fdr>0.1 in NI condition
R_flu2 <- (dplyr::filter(sig_FLU_results, sig_FLU_results$lfsr < .10))
R_NI2 <- (dplyr::filter(NI_results, NI_results$lfsr > .10 ))
Flu_specific <- intersect(R_flu2$feature, R_NI2$feature)
#Number of Flu specific.
print(length(Flu_specific))

## get the rows for the features that are condition specific
NI_specific <- as.vector(NI_specific)
Flu_specific <- as.vector(Flu_specific)

## create function which collects if there is a overlap between the snp and a TF footprint
get_region <- function(vec, id) {
if(length(.footprints <- which(vec >= footprints$start & vec <= footprints$end & id == footprints$chr)[1])) .footprints else NA
}

if(condition=="NI"){
    sig_QTL <- NI_q[NI_q$feature %in% NI_specific, ]
    print("running NI specific")
} else if (condition=="Flu"){
    sig_QTL <- FLU_q[FLU_q$feature %in% Flu_specific, ]
    print("running Flu specific")
}

dim(sig_QTL)[1]
sig_QTL <- select(sig_QTL, snps)
colnames(sig_QTL) <- c("snps")

## obtain snp location for QTL results
sig_QTL$chr <- str_extract(sig_QTL$snps, "[^_]+")
sig_QTL$chr <- paste0("chr", sig_QTL$chr)
sig_QTL$pos <- as.numeric(sub("^[^_]*_([^_]*).*", "\\1", sig_QTL$snps))

for (cluster_num in 1:200){
    print(cluster_num)
    footprints <- footprints_cluster[footprints_cluster$cluster == cluster_num, ]

    for (chr in 1:length(chromosomes)){
        chr_i<- chromosomes[chr]
        chr_i_QTL <- sig_QTL[sig_QTL$chr == chr_i, ]
        if(dim(chr_i_QTL)[1] > 0){
            chr_i_QTL$overlapping_footprint <- as.character(footprints$TF_motif[mapply(get_region, chr_i_QTL$pos, chr_i)])
            if (chr==1) {
                sig_QTL_with_overlap <- chr_i_QTL
            }else {
                sig_QTL_with_overlap <- rbind(sig_QTL_with_overlap, chr_i_QTL)
            }
        }
    }
    sig_QTL_with_overlap$overlapping_footprint[is.na(sig_QTL_with_overlap$overlapping_footprint)] <- "no_overlap"

    ## create df with one column for each TF motif cluster.
    if(cluster_num==1){
        sig_QTL_all_clusters <- sig_QTL_with_overlap
        colnames(sig_QTL_all_clusters)[4] <- cluster_num
    }else {
        sig_QTL_all_clusters <- cbind(sig_QTL_all_clusters, sig_QTL_with_overlap$overlapping_footprint)
        colnames(sig_QTL_all_clusters)[3+cluster_num] <- cluster_num
    }
}

print(dim(sig_QTL_all_clusters)[1]==dim(sig_QTL)[1])

## write lists
saveRDS(sig_QTL_all_clusters,  file = paste0(out_dir,data_type_i, "_sig_", condition, ".cluster.info.rds"))


###############################################################################
## Run same analysis for the best snps of those QTL that are not significant ##
###############################################################################

if(condition=="NI"){
    nonsig_QTL <- NI_q[NI_q$qvalues >= .10, ]
} else if (condition=="Flu"){
    nonsig_QTL <- FLU_q[FLU_q$qvalues >= .10, ]
}

## if data is WGBS which has ~7 million non-significant QTL, randomly subset 500k
if(data_type_i == "WGBS"){
    nonsig_QTL <- nonsig_QTL[sample(rownames(nonsig_QTL), 500000, replace=FALSE), ]
}

## obtain snp location for QTL results
nonsig_QTL$chr <- str_extract(nonsig_QTL$snps, "[^_]+")
nonsig_QTL$chr <- paste0("chr", nonsig_QTL$chr)
nonsig_QTL$pos <- as.numeric(sub("^[^_]*_([^_]*).*", "\\1", nonsig_QTL$snps))

for (cluster_num in 1:200){
    print(cluster_num)
    footprints <- footprints_cluster[footprints_cluster$cluster == cluster_num, ]

    for (chr in 1:length(chromosomes)){
        chr_i<- chromosomes[chr]
        chr_i_QTL <- nonsig_QTL[nonsig_QTL$chr == chr_i, ]
        chr_i_QTL$overlapping_footprint <- as.character(footprints$TF_motif[mapply(get_region, chr_i_QTL$pos, chr_i)])
        if (chr==1) {
            nonsig_QTL_with_overlap <- chr_i_QTL
        }else {
            nonsig_QTL_with_overlap <- rbind(nonsig_QTL_with_overlap, chr_i_QTL)
        }
    }
    nonsig_QTL_with_overlap$overlapping_footprint[is.na(nonsig_QTL_with_overlap$overlapping_footprint)] <- "no_overlap"

    ## create df with one column for each TF motif cluster.
    if(cluster_num==1){
        nonsig_QTL_all_clusters <- nonsig_QTL_with_overlap
        colnames(nonsig_QTL_all_clusters)[10] <- cluster_num
    }else {
        nonsig_QTL_all_clusters <- cbind(nonsig_QTL_all_clusters, nonsig_QTL_with_overlap$overlapping_footprint)
        colnames(nonsig_QTL_all_clusters)[9+cluster_num] <- cluster_num
    }
}

print(dim(nonsig_QTL_all_clusters)[1]==dim(nonsig_QTL)[1])

## write lists
saveRDS(nonsig_QTL_all_clusters,  file = paste0(out_dir,data_type_i, "_NONsig_", condition, ".cluster.info.rds"))
