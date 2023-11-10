## Prepare input data for Torus and get SNP level priors

# Load packages
suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(doParallel))

# Process the command-line arguments
option_list <- list(
  make_option("--phenotype", action="store", default=NA, type='character',
              help="Phenotype name."),
  make_option("--outdir", action="store", default=NA, type='character',
              help="Directory for SuSiE input data."),
  make_option("--path_torus", action="store", default='torus', type='character',
              help="path to vcftools.")
)

opt <- parse_args(OptionParser(option_list=option_list))
phenotype_name <- opt$phenotype
outdir         <- opt$outdir
path_torus     <- opt$path_torus

## Begins here
outdir <- paste0(outdir, "/", phenotype_name)
dir.create(outdir, showWarnings = F, recursive = T)

cat("Prepare Torus input data for", phenotype_name, "...\n")

## Load QTL summary stats and save in MatrixEQTL format
if(!file.exists(paste0("/project2/xinhe/kevinluo/flu_regulatoryQTLs/QTL_results/", phenotype_name, "/full_QTL_sumstats.rds"))){
  QTL_res <- fread(paste0("/project2/xinhe/kevinluo/flu_regulatoryQTLs/QTL_results/", phenotype_name, "/result_original_rsIDs.txt.gz"), drop = 1) %>% as_tibble()
  QTL_res <- QTL_res %>% separate(snps, c("CHROM", "POS", NA, NA), sep = "_", remove = FALSE) %>% unite("snp_pos", CHROM, POS, sep = ":")
  QTL_res$rsID[which(QTL_res$rsID == ".")] <- QTL_res$snp_pos[which(QTL_res$rsID == ".")]
  saveRDS(QTL_res, paste0("/project2/xinhe/kevinluo/flu_regulatoryQTLs/QTL_results/", phenotype_name, "/full_QTL_sumstats.rds"))
}else{
  QTL_res <- readRDS(paste0("/project2/xinhe/kevinluo/flu_regulatoryQTLs/QTL_results/", phenotype_name, "/full_QTL_sumstats.rds"))
}

any(grep("\\.", QTL_res$rsID))
any(grep(":", QTL_res$rsID))

sumstat_file <- paste0(outdir, "/", phenotype_name, ".QTL.sumstat.gz")
fwrite(QTL_res[, c("rsID", "gene", "beta", "statistic", "pvalue")], sumstat_file, sep = "\t")

## Load SNP positions from VCF file and save SNP map
# 'zgrep -v "^##" imputed_ALL.chrs.annotated.v2.renamed.vcf.gz | cut -f1-5 | gzip > snp_pos.txt.gz'
snp_pos_file <- "/project2/xinhe/kevinluo/flu_regulatoryQTLs/genotype/202106_imputed_vcf_v2/snp_pos.txt.gz"
snp_map <- fread(snp_pos_file)
colnames(snp_map)[1] <- "CHR"

snp_map <- snp_map[, c("ID", "CHR", "POS")]
snp_map <- snp_map[snp_map$ID %in% QTL_res$rsID, ]
snp_map_file <- paste0(outdir, "/", phenotype_name, ".snp.map.gz")
fwrite(snp_map, snp_map_file, sep = "\t", col.names = F)

## Load the positions of the molecular phenotype and save gene map
coordinates_dir <- "/project2/xinhe/kevinluo/flu_regulatoryQTLs/positions_files"
gene_map <- read.table(file.path(coordinates_dir, paste0(gsub("_.*", "", phenotype_name), "_positions.txt.gz")), header = T, stringsAsFactors = F)
colnames(gene_map) <- c("gene", "chr", "start", "end")

gene_map <- gene_map[gene_map$gene %in% QTL_res$gene, ]
gene_map$chr <- as.integer(gsub("chr", "", gene_map$chr))
gene_map_file <- paste0(outdir, "/", phenotype_name, ".gene.map.gz")
fwrite(gene_map, gene_map_file, sep = "\t", col.names = F)

# Compute the distance of SNP to nearest gene
QTL_res_tmp <- QTL_res
QTL_res_tmp$snp_pos <- snp_map$POS[match(QTL_res_tmp$rsID, snp_map$ID)]
m_gene <- match(QTL_res_tmp$gene, gene_map$gene)
QTL_res_tmp$gene_start <- gene_map$start[m_gene]
QTL_res_tmp$gene_end <- gene_map$end[m_gene]

QTL_res_tmp$dist <- pmin(abs(QTL_res_tmp$snp_pos - QTL_res_tmp$gene_start),
                           abs(QTL_res_tmp$snp_pos - QTL_res_tmp$gene_end))

QTL_res_tmp <- QTL_res_tmp[order(QTL_res_tmp$dist, decreasing = F),]
QTL_res_tmp <- QTL_res_tmp[!duplicated(QTL_res_tmp$rsID),]

snp_map$dist <- QTL_res_tmp[match(snp_map$ID, QTL_res_tmp$rsID), "dist"]
summary(snp_map$dist)
snp_map$dist[is.na(snp_map$dist)] <- 1e6

## Create an annotation file with distance to peak as annotation (categorical)
annot <- data.frame(SNP = snp_map$ID , dist_d = 0)
annot$dist_d[which(snp_map$dist < 500)] <- 1
annot$dist_d[which(snp_map$dist >= 500 & snp_map$dist < 1000)] <- 2
annot$dist_d[which(snp_map$dist >= 1000 & snp_map$dist < 2000)] <- 3
annot$dist_d[which(snp_map$dist >= 2000 & snp_map$dist < 5000)] <- 4
annot$dist_d[which(snp_map$dist >= 5000 & snp_map$dist < 10000)] <- 5

table(annot$dist_d)

annot_file <- paste0(outdir, "/", phenotype_name, ".dist.annot.gz")
fwrite(annot, annot_file, sep = "\t", col.names = T)

## Run Torus

cat("Run Torus with distance based annotations ...\n")

cat("Estimate enrichment and get SNP level priors ...\n")

tmpdir <- tempdir()
setwd(tmpdir)

prior_dir <- paste0(tmpdir, "/", phenotype_name, ".dist.prior")
system(paste("rm -rf", prior_dir))

torus_cmd <- paste(path_torus,
                   "-d", sumstat_file,
                   "-annot", annot_file,
                   "-dump_prior", prior_dir,
                   "-est >", paste0(outdir, "/", phenotype_name, ".dist.enrichment.est"))

print(torus_cmd)
system(torus_cmd)

# Load and combine priors
registerDoParallel(10)
cat('Using', getDoParWorkers(), 'cores in parallel... \n')

# prior_files <- list.files(prior_dir)
# cat("Load and combine", length(prior_files), "prior results ...\n")
#
# all_priors <- foreach( prior_file=prior_files ) %dopar% {
#   prior <- fread(file.path(prior_dir, prior_file))
#   colnames(prior) <- c("SNP", "prior_prob")
#   prior
# }
# names(all_priors) <- gsub("\\.prior", "", prior_files)
# saveRDS(all_priors, paste0(outdir, "/", phenotype_name, ".dist.priors.rds"))

sig_genes <- unique(QTL_res$gene[QTL_res$FDR < 0.1])

cat("Load and combine", length(sig_genes), "prior files ...\n")

priors <- foreach( iGene=sig_genes ) %dopar% {
  prior <- fread(file.path(prior_dir, paste0(iGene, ".prior")))
  colnames(prior) <- c("SNP", "prior_prob")
  prior
}
names(priors) <- sig_genes
saveRDS(priors, paste0(outdir, "/", phenotype_name, ".dist.priors.rds"))


stopImplicitCluster()


cat("Torus results saved at:", outdir, "\n")
