## Prepare input data for SuSiE QTL fine-mapping

# Load packages
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(doParallel))
suppressMessages(library(susieR))
suppressMessages(library(vcfR))
source("../code/susie_finemapping_utils.R")

# Process the command-line arguments
option_list <- list(
  make_option("--phenotype", action="store", default=NA, type='character',
              help="Phenotype name."),
  make_option("--vcf", action="store", default=NA, type='character',
              help="Genotype VCF file."),
  make_option("--outdir", action="store", default=NA, type='character',
              help="Directory for SuSiE input data."),
  make_option("--ncores", action="store", default=6, type='integer',
              help="Number of CPUs to run in parallel."),
  make_option("--nbatches", action="store", default=5, type='integer',
              help="Number of gene batches to run in parallel."),
  make_option("--ibatch", action="store", default=1, type='integer',
              help="Run for the ith batch of genes."),
  make_option("--path_vcftools", action="store", default='vcftools', type='character',
              help="path to vcftools.")
)

opt <- parse_args(OptionParser(option_list=option_list))
phenotype_name <- opt$phenotype
vcf_file       <- opt$vcf
outdir         <- opt$outdir
num_cores      <- opt$ncores
num_batches    <- opt$nbatches
ibatch         <- opt$ibatch
path_vcftools  <- opt$path_vcftools

# ## Example settings
# phenotype_name <- "RNAseq_NI"
# vcf_file <- "/project2/xinhe/kevinluo/flu_regulatoryQTLs/genotype/202106_imputed_vcf_v2/imputed_ALL.chrs.annotated.v2.renamed.vcf.gz"
# outdir <- "/scratch/midway2/kaixuan/flu_regulatoryQTLs/susie/input"
# num_cores <- 8
# num_batches <- 5
# ibatch <- 1
# path_vcftools <- "~/softwares/vcftools/bin/vcftools"

## Begins here
coordinates_dir <- "/project2/xinhe/kevinluo/flu_regulatoryQTLs/positions_files"
phenotype_dir <- "/project2/xinhe/kevinluo/flu_regulatoryQTLs/expression_matrices"

susie_input_dir <- paste0(outdir, "/", phenotype_name)
dir.create(susie_input_dir, showWarnings = F, recursive = T)

cat("Phenotype:", phenotype_name, "\n")
cat("Genotype:", vcf_file, "\n")
cat("Directory of SuSiE input data:", susie_input_dir, "\n")

## Load QTL data
cat("Loading", phenotype_name, "QTL full individual data ... \n")

## Load QTL summary stats
if(!file.exists(paste0("/project2/xinhe/kevinluo/flu_regulatoryQTLs/QTL_results/", phenotype_name, "/full_QTL_sumstats.rds"))){
  QTL_res <- fread(paste0("/project2/xinhe/kevinluo/flu_regulatoryQTLs/QTL_results/", phenotype_name, "/result_original_rsIDs.txt.gz"), drop = 1) %>% as_tibble()
  QTL_res <- QTL_res %>% separate(snps, c("CHROM", "POS", NA, NA), sep = "_", remove = FALSE) %>% unite("snp_pos", CHROM, POS, sep = ":")
  QTL_res$rsID[which(QTL_res$rsID == ".")] <- QTL_res$snp_pos[which(QTL_res$rsID == ".")]
  saveRDS(QTL_res, paste0("/project2/xinhe/kevinluo/flu_regulatoryQTLs/QTL_results/", phenotype_name, "/full_QTL_sumstats.rds"))
}else{
  QTL_res <- readRDS(paste0("/project2/xinhe/kevinluo/flu_regulatoryQTLs/QTL_results/", phenotype_name, "/full_QTL_sumstats.rds"))
}

## Select genes with FDR < 10%
sig_genes <- unique(QTL_res$gene[which(QTL_res$FDR < 0.1)])
QTL_res <- QTL_res[QTL_res$gene %in% sig_genes, ]
cat("Selected", length(sig_genes), "genes with FDR < 10%. \n")

## Load molecular phenotype data
phenotype_data <- read.table(paste0(phenotype_dir, "/", phenotype_name, ".expression.txt.gz"), sep = "\t", header = T, stringsAsFactors = F)
phenotype_data <- data.frame(gene = rownames(phenotype_data), phenotype_data)
phenotype_data <- phenotype_data %>% filter(gene %in% QTL_res$gene)

## Load coordinates of the molecular phenotype
phenotype_prefix <- gsub("_.*", "", phenotype_name)
coordinates <- read.table(file.path(coordinates_dir, paste0(phenotype_prefix, "_positions.txt.gz")), header = T, stringsAsFactors = F)
colnames(coordinates) <- c("gene", "chr", "start", "end")
phenotype_data <- merge(coordinates, phenotype_data, by = "gene")

phenotype_mat <- t(as.matrix(phenotype_data[, !colnames(phenotype_data) %in% c("gene", "chr", "start", "end")]))
colnames(phenotype_mat) <- phenotype_data$gene

## Covariate matrix
covariate_mat <- NULL # PCs and other covariates were already adjusted.

## create gene ID map
gene_IDs <- sort(unique(QTL_res$gene))

gene_map <- data.frame(gene_name = paste0("Gene",1:length(gene_IDs)), gene_ID = gene_IDs, row.names = gene_IDs, stringsAsFactors = FALSE)
saveRDS(gene_map, file = file.path(susie_input_dir, "gene_map.rds") )

## QTL mapping window size
window_size <- ifelse(grepl("WGBS", phenotype_name), 5000, 1e5)
cat("window size =", window_size, "\n")

## begins preparing the data
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat('Using', getDoParWorkers(), 'cores in parallel... \n')

## Split genes into batches
cat(length(gene_IDs), "genes in total.\n")

split_gene_batches <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

gene_batches <- split_gene_batches(gene_IDs, num_batches)

genes_ibatch <- gene_batches[[ibatch]]

QTL_res <- dplyr::filter(QTL_res, gene %in% genes_ibatch)

cat("Prepare Susie input for", length(genes_ibatch), "genes in batch", ibatch, "...\n")

### check for genes without input data

genes_todo <- foreach( iGene = genes_ibatch, .combine = c )%dopar%{
  if( !file.exists( paste0(susie_input_dir, "/", iGene, ".QTL.full.data.rds") ) ){
    return(iGene)
  }
}

cat(length(genes_todo), "genes to be done. \n")

if(length(genes_todo) > 0){

  ans <- foreach( iGene=genes_todo, .combine = function(a,b){ NULL } ) %dopar% {

    library(tidyverse)

    cat("Prepare Susie input data for: ",iGene, "\n")

    # phenotype vector: n x 1
    y <- phenotype_mat[, iGene]

    gene_bed <- phenotype_data[phenotype_data$gene == iGene, c("chr", "start", "end", "gene")]
    gene_bed$chr <- gsub("chr", "", gene_bed$chr)

    # genotype matrix: n x p, p: number of SNPs
    gene_QTL_res <- QTL_res %>% dplyr::filter(gene == iGene) %>%
      separate(snps, c("CHROM", "POS", NA, NA), sep = "_", remove = FALSE) %>%
      unite("snp_pos", CHROM, POS, sep = ":")

    gene_QTL_res$rsID[which(gene_QTL_res$rsID == ".")] <- gene_QTL_res$snp_pos[which(gene_QTL_res$rsID == ".")]
    snpID_list <- gene_QTL_res$rsID
    cat(length(snpID_list), 'SNPs in this gene window. \n')
    geno.l <- get_genotype_dosage(vcf_file, gene_bed, window_size, snpID_list = snpID_list, path_vcftools = path_vcftools)

    if( is.null(geno.l$geno) )  {
      ## skip this gene if no SNPs were extracted from the VCF file
      cat("No SNPs extracted for this gene. \n")
      saveRDS(NULL, file = paste0(susie_input_dir, "/", iGene, ".QTL.full.data.rds"))
    }else{
      # # extract genotypes using SNP positions
      # gene_QTL_res <- gene_QTL_res %>% separate(snps, c("CHROM", "POS", NA, NA), sep = "_", remove = FALSE)
      # snpPOS_list <- gene_QTL_res[, c("CHROM", "POS")]
      # cat(nrow(snpPOS_list), 'SNPs in this gene window. \n')
      # geno.l <- get_genotype_dosage(vcf_file, gene_bed, window_size, snpPOS_list = snpPOS_list, path_vcftools = path_vcftools)

      X <- t(as.matrix(geno.l$geno, drop= F))
      dim(X)

      ## match the order of individual names
      X <- as.matrix(X[names(y), ])
      colnames(X) <- rownames(geno.l$geno)

      if(!is.null(covariate_mat)){
        covariates <- covariate_mat[names(y), ]
      }else{
        covariates <- NULL
      }

      snps_withNAs <- apply(X, 2, function(x) anyNA(x))
      if(any(snps_withNAs)){
        cat('Remove', length(snps_withNAs), 'variants with missing values. \n')
        X <- X[, !snps_withNAs]
      }

      ## select common SNPs
      common_SNPs <- intersect(colnames(X), gene_QTL_res$rsID)
      X <- X[, match(common_SNPs, colnames(X)), drop = F]
      gene_QTL_res <- gene_QTL_res[match(common_SNPs, gene_QTL_res$rsID), ]

      cat("Dimension of X matrix: ", nrow(X), "x", ncol(X), "\n")

      cat("length of y: ", length(y), "\n")

      cat("Dimension of covariates matrix: \n")
      print(dim(covariates))

      saveRDS(list(X = X, y = y, covariates = covariates),
              file = paste0(susie_input_dir, "/", iGene, ".QTL.full.data.rds"))

    }
  }

}


genes_todo <- foreach( iGene = genes_ibatch, .combine = c )%dopar%{
  if( !file.exists( paste0(susie_input_dir, "/", iGene, ".QTL.full.data.rds") ) ){
    return(iGene)
  }
}

cat(length(genes_todo), "genes without input data. \n")

# Stop cluster
stopCluster(cl)
