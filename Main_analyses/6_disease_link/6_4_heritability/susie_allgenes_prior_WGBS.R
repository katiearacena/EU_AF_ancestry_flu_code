## Fine-mapping QTLs using SuSiE
setwd("~/projects/flu_regulatoryQTLs/code")
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
  make_option("--L", action="store", default=10, type='integer',
              help="Number of causal effects per locus (default = 10)."),
  make_option("--susie_prior", action="store", default="uniform", type='character',
              help="Type of prior for susie (dist or uniform)."),
  make_option("--torus_dir", action="store", default=NA, type='character',
              help="Directory for Torus priors."),
  make_option("--susie_input_dir", action="store", default=NA, type='character',
              help="Directory for SuSiE input."),
  make_option("--outdir", action="store", default=NA, type='character',
              help="Directory for SuSiE output."),
  make_option("--ncores", action="store", default=6, type='integer',
              help="Number of CPUs to run in parallel.")
)

opt <- parse_args(OptionParser(option_list=option_list))
phenotype_name  <- opt$phenotype
L               <- opt$L
susie_prior     <- opt$susie_prior
torus_dir       <- opt$torus_dir
susie_input_dir <- opt$susie_input_dir
outdir          <- opt$outdir
num_cores       <- opt$ncores

# # Example settings
# phenotype_name <- "WGBS_Flu_full_E3"
# L <- 3
# susie_prior <- "dist"
# torus_dir <- "/scratch/midway2/kaixuan/flu_regulatoryQTLs/torus/WGBS_Flu_full_E7"
# susie_input_dir <- "/scratch/midway2/kaixuan/flu_regulatoryQTLs/susie/input/WGBS_Flu_full_E7"
# outdir <- "/project2/xinhe/kevinluo/flu_regulatoryQTLs/susie/output/WGBS_Flu_full_E7"
# num_cores <- 6

## begins here
registerDoParallel(num_cores)
cat('Using', getDoParWorkers(), 'cores in parallel... \n')

## load QTL results
susie_output_dir <- paste0(outdir, "/susie_", susie_prior, "_prior")
dir.create(susie_output_dir, showWarnings = F, recursive = T)


cat("Phenotype name:", phenotype_name, "\n")
cat("SuSiE L=", L, "\n")
cat("SuSiE input data dir:", susie_input_dir, "\n")
cat("SuSiE output dir:", susie_output_dir, "\n")

## Load susie input data
cat("Loading", phenotype_name, "input data ... \n")

gene_map <- readRDS(file.path(susie_input_dir, "gene_map.rds"))

### check input data
gene_IDs <- gene_map$gene_ID
cat(length(gene_IDs), "genes in the gene map. \n")

genes_without_input <- foreach( iGene = gene_IDs, .combine = c )%dopar%{
  if( !file.exists( paste0(susie_input_dir, "/", iGene, ".QTL.full.data.rds") ) ){
    return(iGene)
  }
}

if(length(genes_without_input) > 0){
  cat("Warning:", length(genes_without_input), "genes without input data! \n")
}

genes_to_finemap <- setdiff(gene_IDs, genes_without_input)
cat(length(genes_to_finemap), "genes for fine-mapping. \n")

## Run susie finemapping for all genes

cat("Finemapping: ", length(genes_to_finemap), "genes \n")

cat("Use", susie_prior, "prior ...\n")
if(susie_prior != "uniform"){
  cat("torus dir:", torus_dir, "\n")
  torus_priors <- readRDS(paste0(torus_dir, "/", phenotype_name, ".dist.priors.rds"))
}

if(length(genes_to_finemap) > 0){

  susie_allgenes_res <- foreach( iGene=genes_to_finemap ) %dopar% {

    cat("Susie finemapping for: ",iGene, "\n")

    data <- readRDS(paste0(susie_input_dir, "/", iGene, ".QTL.full.data.rds"))
    if(is.null(data)){
      cat("No input data for", iGene, "\n")
      return(NULL)
    }else{

      if(susie_prior == "dist"){
        ## load prior prob. obtained from Torus
        prior <- torus_priors[[iGene]]
        prior <- prior[match(colnames(data$X), prior$SNP),] ## match SNP order
        prior[is.na(prior$prior_prob)] <- 0

        # Susie with dist prior
        susie_gene_res <- fit_susie_idv(data$X, data$y, data$covariates, prior = prior$prior_prob, L)

      }else{
        # Susie with uniform prior
        susie_gene_res <- fit_susie_idv(data$X, data$y, data$covariates, prior = NULL, L)
      }

      return(susie_gene_res)
    }

  }

}

stopImplicitCluster()

if(length(genes_to_finemap) != length(susie_allgenes_res)){
  stop("Check: unequal number of genes! \n")
}

names(susie_allgenes_res) <- as.character(genes_to_finemap)

saveRDS(susie_allgenes_res, file = paste0(susie_output_dir, "/susie_allgenes_res.rds") )

cat("Finemapped",length(susie_allgenes_res), "genes. \n")

## PVE for all genes
PVE_all <- unlist(sapply(1:length(susie_allgenes_res), function(i) susie_allgenes_res[[i]]$PVE))
cat("PVE for", phenotype_name, "in all genes: \n")
print(summary(PVE_all))


## Combine summary statistics for all gene-SNP pairs
SusieResult <- combine_susie_results(susie_allgenes_res, num_cores = num_cores)
SusieResult <- SusieResult %>% as_tibble() %>% dplyr::rename(rsID = SNP, gene = GENE)

cat("Credible sets:\n")
table(SusieResult$CS)
cat(length(unique(SusieResult$gene[SusieResult$CS != 0])), "genes have credible sets \n")

## Add SNP positions and other information

## Load QTL summary stats
if(!file.exists(paste0("/project2/xinhe/kevinluo/flu_regulatoryQTLs/QTL_results/subset_WGBS/", phenotype_name, "_QTL_sumstats.rds"))){
  QTL_res <- fread(paste0("/project2/xinhe/kevinluo/flu_regulatoryQTLs/QTL_results/subset_WGBS/", phenotype_name, ".txt.gz")) %>% as_tibble()
  QTL_res <- QTL_res %>% separate(snps, c("CHROM", "POS", NA, NA), sep = "_", remove = FALSE) %>% unite("snp_pos", CHROM, POS, sep = ":")
  QTL_res$rsID[which(QTL_res$rsID == ".")] <- QTL_res$snp_pos[which(QTL_res$rsID == ".")]
  saveRDS(QTL_res, paste0("/project2/xinhe/kevinluo/flu_regulatoryQTLs/QTL_results/subset_WGBS/", phenotype_name, "_QTL_sumstats.rds"))
}else{
  QTL_res <- readRDS(paste0("/project2/xinhe/kevinluo/flu_regulatoryQTLs/QTL_results/subset_WGBS/", phenotype_name, "_QTL_sumstats.rds"))
}

# Add other information from the QTL results
SusieResult <- SusieResult %>% left_join(QTL_res, by = c("rsID", "gene")) %>%
  separate(snps, c("chr", "pos", NA, NA), sep = "_", remove = FALSE) %>%
  mutate(pos = as.integer(pos)) %>%
  select(rsID, snps, chr, pos, gene, beta, statistic, pvalue, FDR, PIP, CS)

saveRDS(SusieResult, file = paste0(susie_output_dir, "/susie_", phenotype_name, "_", susie_prior, "_prior_sumstats.rds"))

cat("Done SuSiE finemapping.\n")
