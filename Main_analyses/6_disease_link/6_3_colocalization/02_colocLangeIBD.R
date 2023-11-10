###############################################################################
## Script purpose: COLOC with macrophage QTLs
## Author: Zepeng Mu
## Date: Thu Sep 24 19:03:53 2020
###############################################################################
library(coloc)
library(tidyverse)
library(data.table)
"%&%" <- function(a, b) paste0(a, b)

dataDir <- "/project2/lbarreiro/users/katie/EU_AF_ancestry_flu/admixture/final_outputs"
resultsDir <- "/project2/lbarreiro/users/katie/EU_AF_ancestry_flu/paper_pipeline/Outputs/3_QTL_mapping/SNP-QTL_mapping"

## Load GWAS full summary statistics
## Format Chr, SNP and BP columns properly to match QTL data
gwas <- fread("/project2/yangili1/zpmu/GWAS_loci/ibd_build37_59957_20161107.txt.gz")
gwas <- gwas %>%
  mutate(SNP = sapply(strsplit(gwas$MarkerName, "_"), "[[", 1),
         CHR = sapply(strsplit(gwas$MarkerName, ":"), "[[", 1),
         BP = sapply(strsplit(SNP, ":"), "[[", 2))

gwas$CHR <- as.integer(gwas$CHR)
gwas$BP <- as.integer(gwas$BP)

## Load GWAS loci table generated from previous script
## Remove HLA region
loci <- read.table("/project2/yangili1/zpmu/GWAS_loci/lange_ibd_loci_1e-5.txt",
                   header = T) %>%
  filter(!(Chr == 6 & Position > 25e6 & Position < 35e6))

## Iterate through each phenotype and condition
conditions <- c("Flu", "NI")
pheno <- c("RNAseq", "ATACseq", "H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "WGBS")
pheno <- "WGBS"

for (ana in pheno) {
  ## Create output dataframe
  outDf <- data.frame(condition = character(), gene = character(), colocSnp = character(),
                      lociName = character(), hit2 = character(), hit1 = character(),
                      nsnps = double(), PP.H0.abf = double(), PP.H1.abf = double(),
                      PP.H2.abf = double(), PP.H3.abf = double(), PP.H4.abf = double(),
                      best1 = character(), best2 = character(), best4 = character(),
                      hit1.margz = double(), hit2.margz = double())
  
  ## Smaller window size used for WGBS
  if (ana == "WGBS") {
    windowSize <- 5000
  } else {
    windowSize <- 100e3
  }
  
  for (cond in conditions) {
    ## Load permutation QTL results and extract significant QTLs
    permQtl.name <- "results_best_SNPs_with_qval.txt"
    permQtl <- fread(file.path(resultsDir, ana, cond, ana%&%"_"%&%cond, permQtl.name)) %>%
      filter(qvalues < 0.1) %>%
      separate(snps, into = c("Chr", "Position", "a1", "a2"), sep = "_", remove = F) %>%
      mutate(sid = Chr%&%":"%&%Position,
             Chr = as.integer(Chr),
             Position = as.integer(Position)) %>%
      select(pid = gene, sid, Chr, Position, qvalues, slope = beta, snps)
    
    ## Load nominal QTL results
    nomQtl.name <- "result_original_rsIDs.txt"
    nomQtl <- fread(file.path(resultsDir, ana, cond, ana%&%"_"%&%cond, "raw_results", nomQtl.name)) %>%
      separate(snps, into = c("Chr", "Position", "a1", "a2"), sep = "_", remove = F) %>%
      mutate(sid = Chr%&%":"%&%Position,
             Chr = as.integer(Chr),
             Position = as.integer(Position)) %>%
      select(pid = gene, sid, Chr, Position, npval = pvalue, slope = beta, snps)
    
    ## Load MAF data in QTL study
    qtlMaf <- fread(
      file.path(dataDir, "bi-allelic-only.35Samples.noIPSC.hc.vqsr.filtered.rs.ID.annotated.v4.frq"),
      col.names = c("CHR", "POS", "N_All", "N_CHR", "freq1", "freq2"),
      skip = 1,
      header = F
    ) %>%
      separate(freq1, c("a1", "freqa1"), sep = ":") %>%
      separate(freq2, c("a2", "freqa2"), sep = ":") %>%
      mutate(SNP = CHR%&%":"%&%POS,
             freqa1 = as.double(freqa1),
             freqa2 = as.double(freqa2),
             MAF = case_when(freqa1 < freqa2 ~ freqa1, T ~ freqa2)) %>%
      select(CHR, SNP, MAF)
    
    ## Iterate through each chromosome
    for (i in unique(loci$Chr)) {
      chrLoci <- loci %>% filter(Chr == i)
      chrPermQtl <- permQtl %>% filter(Chr == i)
      chrGwas <- gwas %>% filter(CHR == i)
      chrNomQtl <- nomQtl %>% filter(Chr == i)
      chrMafQtl <- qtlMaf %>% filter(CHR == i)
      
      ## Iterate through each locus
      for (j in 1:nrow(chrLoci)) {
        lociName <- chrLoci$SNPID[j]
        lociPos <- chrLoci$Position[j]
        upStream <- lociPos + windowSize
        dnStream <- lociPos - windowSize
        tmpGene <- filter(chrPermQtl, Position < upStream & Position > dnStream)$pid
        
        if (length(tmpGene) > 0) {
          for (g in tmpGene) {
            tmpQtl <- chrNomQtl %>%
              filter(pid == g & Position < upStream & Position > dnStream & !is.na(npval))
            
            ## Format GWAS and QTL for COLOC input
            tmpQtl <- left_join(tmpQtl, select(chrMafQtl, SNP, MAF), by = c("sid" = "SNP"))
            tmpQtl$abs.Z <- qnorm(tmpQtl$npval / 2, lower.tail = F)
            tmpQtl$var <- (tmpQtl$slope / tmpQtl$abs.Z) ^ 2
            tmpQtl <- tmpQtl %>% filter(!is.na(var) & var > 0)
            tmpGwas <- chrGwas %>% filter(BP < upStream & BP > dnStream)
            
            overlap.snp <- intersect(tmpQtl$sid, tmpGwas$SNP)
            tmpQtl <- tmpQtl %>% filter(sid %in% overlap.snp)
            tmpGwas <- tmpGwas %>% filter(SNP %in% overlap.snp)
            
            ## Remove duplicated SNP IDs
            dup.qtl <- which(duplicated(tmpQtl$sid))
            if (length(dup.qtl) > 0) {
              tmpQtl <- tmpQtl[-dup.qtl, ]
            }
            
            dup.gwas <- which(duplicated(tmpGwas$SNP))
            if (length(dup.gwas) > 0) {
              tmpGwas <- tmpGwas[-dup.gwas, ]
            }
            
            overlap.snp <- intersect(tmpQtl$sid, tmpGwas$SNP)
            tmpQtl <- tmpQtl %>% filter(sid %in% overlap.snp)
            tmpGwas <- tmpGwas %>% filter(SNP %in% overlap.snp)
            
            if (nrow(tmpGwas) > 1) {
              print(g)
              
              my.res <- coloc.signals(
                p12 = 1e-5,
                dataset1 = list(
                  N = 36,
                  beta = tmpQtl$slope,
                  varbeta = tmpQtl$var,
                  type = "quant",
                  MAF = tmpQtl$MAF,
                  snp = tmpQtl$sid
                ),
                dataset2 = list(
                  beta = tmpGwas$Effect,
                  varbeta = tmpGwas$StdErr ^ 2,
                  type = "cc",
                  s = 0.417666,
                  N = 25042 + 34915,
                  snp = tmpGwas$SNP
                )
              )
              
              ## Get COLOC SNP (i.e. SNP with largest PP4)
              colocSnp <- my.res$results$snp[which.max(my.res$results$SNP.PP.H4)]
              for (n in seq_len(nrow(my.res$summary))) {
                outDf <- rbindlist(list(outDf, cbind(cond, g, colocSnp, lociName, my.res$summary[n, ])),
                                   use.names = F)
              }
            }
          }
        }
      }
    }
  }
  
  rm(nomQTL)
  out.dir <- "/project2/yangili1/zpmu/barreiroLab/katie/results/coloc"
  
  fwrite(
    outDf,
    file.path(out.dir, ana%&%"_lange_ibd_noHLA.txt"),
    quote = F,
    row.names = F,
    sep = "\t"
  )
}
