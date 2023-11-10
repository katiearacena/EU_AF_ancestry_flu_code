###############################################################################
## Script purpose: Identify GWAS loci with 1Mb size
## Author: Zepeng Mu
###############################################################################
library(tidyverse)
library(data.table)
"%&%" <- function(a, b) paste0(a, b)
size <- 5e5 ## Half GWAS locus size


## Load GWAS
gwas <- fread("/project2/yangili1/zpmu/GWAS_loci/ibd_build37_59957_20161107.txt.gz")
gwas <- gwas %>%
  mutate(SNPID = sapply(strsplit(gwas$MarkerName, "_"), "[[", 1),
         Chr = sapply(strsplit(gwas$MarkerName, ":"), "[[", 1),
         Position = sapply(strsplit(SNPID, ":"), "[[", 2)) %>%
  select(SNPID, Chr, Position, Pval = P.value)

gwas$Chr <- as.integer(gwas$Chr)
gwas$Position <- as.integer(gwas$Position)
gwas$Pval <- as.double(gwas$Pval)

outDf <- data.frame(SNPID = NA, Chr = NA, Position = NA, Pval = NA)

for (i in 1:22) {
  tmpGwas <- gwas %>%
    filter(Chr == i & Pval < 1e-5) %>%
    arrange(Pval)
  while (nrow(tmpGwas) > 1) {
    outDf <- rbind(outDf, tmpGwas[1, ])
    tmpPos <- tmpGwas$Position[1]
    tmpGwas <- filter(tmpGwas, Position > tmpPos + size | Position < tmpPos - size)
  }
}

outDf <- outDf[complete.cases(outDf), ]

fwrite(outDf, "/project2/yangili1/zpmu/GWAS_loci/lange_ibd_loci_1e-5.txt",
            quote = F, sep = "\t", row.names = F)
