argv <- commandArgs(trailingOnly = TRUE)
source("../../scripts/GTEx_Tissue_Wide_CV_elasticNet.R")

"%&%" <- function(a,b) paste(a, b, sep = "")

tis <- argv[1]
chrom <- argv[2]
alpha <- as.numeric(argv[3])
window <- as.numeric(argv[4])
condition <- argv[5]


data_dir <- "../../data/intermediate/"

expression_RDS <- data_dir %&% "expression_phenotypes/" %&% tis %&% "_batch.age.corrected.RDS"
geno_file <- data_dir %&% "genotypes/" %&% "SNP_genotypes.chr" %&% chrom %&% ".txt"
gene_annot_RDS <- data_dir %&% "annotations/gene_annotation/" %&% tis %&% "_anno.parsed.RDS"
snp_annot_RDS <- data_dir %&% "annotations/snp_annotation/SNPs.annot.chr" %&% chrom %&% ".RDS"
n_k_folds <- 10
out_dir <- data_dir %&% "model_by_chr/"
snpset <- "ME_snps"

TW_CV_model(expression_RDS, geno_file, gene_annot_RDS, snp_annot_RDS,
    n_k_folds, alpha, out_dir, tis, chrom, snpset, window, condition)
