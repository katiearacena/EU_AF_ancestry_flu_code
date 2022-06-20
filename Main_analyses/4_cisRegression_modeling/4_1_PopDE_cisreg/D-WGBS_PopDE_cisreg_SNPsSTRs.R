########################################
### Load libraries and set functions ###
########################################

library(bsseq)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(qvalue)
library(DSS)
library(mashr)

data <- "WGBS"
# Decide whether to permute data (FALSE) or load previous permutations(TRUE)
reproduce <- "FALSE"
# Perform 10 permutations to use for FDR correction.
iterations <- 5

## Set directory structure
folder = "4_cisregression_modeling"
## Set directory 3 steps above script.
setwd('../../../')

## Create directory structure to save outputs.
system(paste0("mkdir -p Outputs/",folder,"/cisregSNPsSTRs/", data,"/"))
out_dir <- paste0("Outputs/",folder,"/cisregSNPsSTRs/", data,"/")

## Set functions
# This function mean centers the data.
mean_center_clean=function(cols_local,column){
  index=which(colnames(cols_local)==column)
  cols_local[,index]=cols_local[,index]-mean(cols_local[,index])
  return(cols_local)
}
# This function applies mean centers function to desired metadata column. Here we apply it to age.
pretty_up_cols=function(cols){

  cols=refactorize(cols)
  cols=mean_center_clean(cols,"Age")
  return(cols)
}
# This function makes sure that desired metadata columns are factors.
refactorize=function(cols_local){

  cols_local$Genotyping_ID=factor(cols_local$Genotyping_ID,levels=unique(as.character(cols_local$Genotyping_ID)))
  cols_local$Sample=factor(cols_local$Sample)
  cols_local$Condition=factor(cols_local$Condition)
  cols_local$Batch=factor(cols_local$Batch)
  return(cols_local)
}

##############################
### Load data and metadata ###
##############################

#Load metadata
cols = read.table(paste0("Inputs/metadata/", data, "_metadata.txt"))
## Load unsmoothed BSseq object
load("Inputs/counts_matrices/WGBS_filtered.counts.BSobj.RData")

## Load and subset SNP genotypes
## read in SNP genotypes (used in matrixeqtl)
SNPgenotypes = read.table(paste0("Inputs/QTL_mapping/SNP_genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)
SNPgenotypes$snp <- rownames(SNPgenotypes)

## read in top SNPs for NI condition
NI_topSNPs <- read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/NI/", data,"_NI/results_best_SNPs_with_qval.txt"))
## only keep snp and gene columns
NI_topSNPs <- NI_topSNPs[,c("snps", "gene")]
## subset genotypes to pertinent ones
NI_SNPgenotypes_subset <- right_join(SNPgenotypes, NI_topSNPs, by = c("snp" = "snps"))
rownames(NI_SNPgenotypes_subset) <- NI_SNPgenotypes_subset$gene; NI_SNPgenotypes_subset$gene <- NULL; NI_SNPgenotypes_subset$snp <- NULL
## Change missing data from -9 to NA.
NI_SNPgenotypes_subset[NI_SNPgenotypes_subset == -9] <- NA

## prepare for adding genotype to meta data
NI_gene_SNPgenotype <- NI_SNPgenotypes_subset[,colnames(NI_SNPgenotypes_subset)]
## add NI to col name
colnames(NI_gene_SNPgenotype) <- paste0(colnames(NI_gene_SNPgenotype), "_NI")

## read in top SNPs for FLU condition
FLU_topSNPs <- read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/Flu/", data,"_Flu/results_best_SNPs_with_qval.txt"))
## only keep snp and gene columns
FLU_topSNPs <- FLU_topSNPs[, c("snps", "gene")]
## subset genotypes to pertinent ones
FLU_SNPgenotypes_subset <- right_join(SNPgenotypes, FLU_topSNPs,  by = c("snp" = "snps"))
rownames(FLU_SNPgenotypes_subset) <- FLU_SNPgenotypes_subset$gene; FLU_SNPgenotypes_subset$gene <- NULL; FLU_SNPgenotypes_subset$snp <- NULL
## Change missing data from -9 to NA.
FLU_SNPgenotypes_subset[FLU_SNPgenotypes_subset == -9] <- NA

## prepare for adding genotype to meta data
FLU_gene_SNPgenotype <- FLU_SNPgenotypes_subset[,colnames(FLU_SNPgenotypes_subset)]
## add FLU to row name
colnames(FLU_gene_SNPgenotype) <- paste0(colnames(FLU_gene_SNPgenotype), "_Flu")

## merge the genotypes of the best SNP in NI and Flu condition
NI_gene_SNPgenotype$feature <- rownames(NI_gene_SNPgenotype)
FLU_gene_SNPgenotype$feature <- rownames(FLU_gene_SNPgenotype)

both_gene_SNPgenotype <-full_join(NI_gene_SNPgenotype, FLU_gene_SNPgenotype, by="feature")
rownames(both_gene_SNPgenotype) <- both_gene_SNPgenotype$feature; both_gene_SNPgenotype$feature <- NULL


## Load and subset STR genotypes
## read in STR genotypes (used in matrixeqtl)
STRgenotypes = read.table(paste0("Inputs/QTL_mapping/STR_genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)
STRgenotypes$str <- rownames(STRgenotypes)

## read in top STRs for NI condition
NI_topSTRs <- read.table(paste0("Outputs/3_QTL_mapping/STR-QTL_mapping/",data,"/NI/", data,"_NI/results_best_STRs_with_qval.txt"))
## only keep str and gene columns
NI_topSTRs <- NI_topSTRs[, c("snps", "gene")]
## subset genotypes to pertinent ones
NI_STRgenotypes_subset <- right_join(STRgenotypes, NI_topSTRs, by = c("str" = "snps"))
rownames(NI_STRgenotypes_subset) <- NI_STRgenotypes_subset$gene; NI_STRgenotypes_subset$gene <- NULL; NI_STRgenotypes_subset$str <- NULL
## Change missing data from -9 to NA.
NI_STRgenotypes_subset[NI_STRgenotypes_subset == -9] <- NA

## prepare for adding genotype to meta data
NI_gene_STRgenotype <- NI_STRgenotypes_subset[,colnames(NI_STRgenotypes_subset)]
## add NI to col name
colnames(NI_gene_STRgenotype) <- paste0(colnames(NI_gene_STRgenotype), "_NI")

## read in top STRs for FLU condition
FLU_topSTRs <- read.table(paste0("Outputs/3_QTL_mapping/STR-QTL_mapping/",data,"/Flu/", data,"_Flu/results_best_STRs_with_qval.txt"))
## only keep snp and gene columns
FLU_topSTRs <- FLU_topSTRs[, c("snps", "gene")]
## subset genotypes to pertinent ones
FLU_STRgenotypes_subset <- right_join(STRgenotypes, FLU_topSTRs,  by = c("str" = "snps"))
rownames(FLU_STRgenotypes_subset) <- FLU_STRgenotypes_subset$gene; FLU_STRgenotypes_subset$gene <- NULL; FLU_STRgenotypes_subset$str <- NULL
## Change missing data from -9 to NA.
FLU_STRgenotypes_subset[FLU_STRgenotypes_subset == -9] <- NA

## prepare for adding genotype to meta data
FLU_gene_STRgenotype <- FLU_STRgenotypes_subset[,colnames(FLU_STRgenotypes_subset)]
## add FLU to col name
colnames(FLU_gene_STRgenotype) <- paste0(colnames(FLU_gene_STRgenotype), "_Flu")

## merge the genotypes of the best STR in NI and Flu condition
NI_gene_STRgenotype$feature <- rownames(NI_gene_STRgenotype)
FLU_gene_STRgenotype$feature <- rownames(FLU_gene_STRgenotype)

both_gene_STRgenotype <-full_join(NI_gene_STRgenotype, FLU_gene_STRgenotype, by="feature")
rownames(both_gene_STRgenotype) <- both_gene_STRgenotype$feature; both_gene_STRgenotype$feature <- NULL


######################################
## Format reads matrices & metadata ##
######################################

#subset to only include those samples for which there is a 1 in the "PopDEset" column 
cols=cols[which(cols$PopDE_set==1),]
#ensure that the Sample_ID are the rownames of the cols dataframe
rownames(cols)=cols$Sample
#ensure that the cols dataframe is ordered alphabetically (based on sample_id)
cols=cols[order(rownames(cols)),]
## Mean center age & refactorize 
cols=pretty_up_cols(cols)
## Set NI as first level so that it is used as the reference.
cols$Condition = factor(cols$Condition, levels=c("NI","Flu"))

## Check to make sure the number of samples in metadata and WGBS reads data matches.
dim(cols)[1] == dim(BSobj.fit.nolow)[2]
dim(BSobj.fit.nolow)[1]==7463164 # number of cpg sites there should be


## Filter based on coverage
#To minimize noise in methylation estimates due to low-coverage data, we restricted analysis to CpG sites
# with coverage of ≥4 sequence reads in at least half of the samples in each condition. DSS accounts for low coverage
# in the model so it is not necessary before.
BS.cov <- getCoverage(BSobj.fit.nolow)
keepLoci.ex <- which(rowSums(BS.cov[, BSobj.fit.nolow$Condition == "Flu"] >= 4) > 17 &
                       rowSums(BS.cov[, BSobj.fit.nolow$Condition == "NI"] >= 4) > 17 )
length(keepLoci.ex)
BSobj.fit.nolow <- BSobj.fit.nolow[keepLoci.ex,]

## remove sites with no variation
reads_whole <-  getMeth(BSobj.fit.nolow, type="raw")
no.var <- rownames(reads_whole[rowVars(reads_whole, na.rm=TRUE) == 0, ])
reads_whole.lim <- subset(reads_whole, !(rownames(reads_whole) %in% no.var))
BSobj.fit.nolow.var <- BSobj.fit.nolow[rownames(reads_whole.lim),]

# remove X, Y, MT and contig features from the matrix.
positions <- read.table(paste0("Inputs/QTL_mapping/", data, "_positions.txt"))
autosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
"chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
positions_filtered <- positions[positions$chromosome %in% autosomes,]
BS.fit <- BSobj.fit.nolow.var[rownames(BSobj.fit.nolow.var) %in% positions_filtered$Gene_ID, ]


##################################################
### Model DE within each condition (real data) ###
##################################################

## subset on features that were originally popDE
## load those popde results
orig_popDE_NI <- read.table(paste0("Outputs/2_ancestry_effects_modeling/PopDE_analysis/", data, "/results/popDE_NI.txt"))
orig_popDE_FLU <- read.table(paste0("Outputs/2_ancestry_effects_modeling/PopDE_analysis/", data, "/results/popDE_FLU.txt"))

## get sig.
sig_NI <- orig_popDE_NI[orig_popDE_NI$qvals_NI < .10, ]
sig_Flu <- orig_popDE_FLU[orig_popDE_FLU$qvals_Flu < .10, ]

## get features that are sig. in either
sig_in_either <- c(rownames(sig_NI), rownames(sig_Flu))
sig_in_either <- unique(sig_in_either)
length(sig_in_either)

## filter down bsseq object that are sig in either
BS.fit <- BS.fit[rownames(BS.fit) %in% sig_in_either, ]
colnames(BS.fit) <- sub("_dup", "", colnames(BS.fit))

##################################
#### Add genotype to metadata  ###
##################################

both_variants_genes <- intersect(rownames(FLU_STRgenotypes_subset), rownames(FLU_SNPgenotypes_subset))
genes <- both_variants_genes[which(both_variants_genes %in% sig_in_either)]

for(j in 1:length(genes)){
  tryCatch({
    if(j%%10 == 0) print(j)
    gene <- genes[j]


    meta_data_loop <- cols
    rownames(meta_data_loop) <- sub("_dup", "", rownames(meta_data_loop))
    meta_data_loop$indiv_ID <- rownames(meta_data_loop)
  
  
    ## add genotype into meta data
    meta_data_loop$SNP <- as.numeric(both_gene_SNPgenotype[gene, rownames(meta_data_loop)])
    meta_data_loop$STR <- as.numeric(both_gene_STRgenotype[gene, rownames(meta_data_loop)])

    ###############
    ## WITH DSS  ##
    ###############

    ## Set design matrix
    design <- select(meta_data_loop, Condition, Admixture, Age, Batch, SNP, STR)
    ## remove inds if there's an NA in SNP genotype column
    design <- design[!is.na(design$SNP), ]
    ## remove inds if there's an NA in STR genotype column
    design <- design[!is.na(design$STR), ]
    ## remove dup from BSseq object column  names 
    BS2 <- BS.fit[ , colnames(BS.fit) %in% rownames(design)]

    ## run the linear model using DSS
    lm = DMLfit.multiFactor(BS2, design=design, formula=~Condition + Condition:Admixture + Condition:SNP + Condition:STR + Condition:Age + Condition:Batch)

    ## extract the admixture coefficients
    NI_results = DMLtest.multiFactor(lm, coef="ConditionNI:Admixture")
    FLU_results = DMLtest.multiFactor(lm, coef="ConditionFlu:Admixture")

    rownames(NI_results) <- paste0(NI_results$chr, "_", NI_results$pos)
    rownames(FLU_results) <- paste0(FLU_results$chr, "_", FLU_results$pos)

    ## for each gene, store each output separately
    results_NI=NI_results[rownames(NI_results) %in% gene , ]
    results_FLU=FLU_results[rownames(FLU_results) %in% gene , ]

    if(j == 1){
      NI_outs <- results_NI
      FLU_outs <- results_FLU
      write.table(NI_outs, paste0(out_dir, "popDE_NI.real.txt"))
      write.table(FLU_outs, paste0(out_dir, "popDE_FLU.real.txt"))
      
    }else{
      NI_outs <- rbind(NI_outs, results_NI)
      FLU_outs <- rbind(FLU_outs, results_FLU)
      write.table(NI_outs, paste0(out_dir, "popDE_NI.real.txt"))
      write.table(FLU_outs, paste0(out_dir, "popDE_FLU.real.txt"))
    }
    
  }, error=function(e){
    #err_genes = c(err_genes, gene)
    #cat("ERROR :",conditionMessage(e), "\n")
  })
}

  
#results_NI <- NI_outs
#results_FLU <- FLU_outs

#write.table(results_NI, paste0(out_dir, "popDE_NI.real.txt"))
#write.table(results_FLU, paste0(out_dir, "popDE_FLU.real.txt"))


###########################################################
### Do permutation tests to define null model for PopDE ###
###########################################################

if(reproduce)
{
  ### to reproduce results previously generated load the permuted data.
  shuffled_pvals_adm_NI=read.table(paste0(out_dir, "/permuted_pvalues/popDE_NI.txt"))
  shuffled_pvals_adm_FLU=read.table(paste0(out_dir, "/permuted_pvalues/popDE_FLU.txt"))
  ##if reproduce is set to TRUE this part below will be skipped. If reproduce = FALSE then you the code below will permute the data.
}else{

  design$SNP <- NULL
  design$STR <- NULL
  design_random=design

  for(iter in 1:iterations)
  {
    if(iter%%1==0)print(iter)
    design_random$Admixture=sample(design_random$Admixture)

    for(j in 1:length(genes)){
      tryCatch({
        if(j%%10 == 0) print(j)
        gene <- genes[j]

        meta_data_loop <- design_random
        rownames(meta_data_loop) <- sub("_dup", "", rownames(meta_data_loop))
        meta_data_loop$indiv_ID <- rownames(meta_data_loop)
      
      
        ## add genotype into meta data
        meta_data_loop$SNP <- as.numeric(both_gene_SNPgenotype[gene, rownames(meta_data_loop)])
        meta_data_loop$STR <- as.numeric(both_gene_STRgenotype[gene, rownames(meta_data_loop)])

        ###############
        ## WITH DSS  ##
        ###############

        ## Set design matrix
        design_perm <- select(meta_data_loop, Condition, Admixture, Age, Batch, SNP, STR)
        ## remove inds if there's an NA in genotype column
        design_perm <- design_perm[!is.na(design_perm$genotype), ]
        ## remove dup from BSseq object column  names 
        BS2 <- BS.fit[ , colnames(BS.fit) %in% rownames(design_perm)]

        lm_rand = DMLfit.multiFactor(BS2, design=design_perm, formula=~Condition + Condition:Admixture + Condition:SNP + Condition:STR + Condition:Age + Condition:Batch)
      
        ## extract the admixture coefficients
        NI_perm = DMLtest.multiFactor(lm_rand, coef="ConditionNI:Admixture")
        FLU_perm = DMLtest.multiFactor(lm_rand, coef="ConditionFlu:Admixture")

        rownames(NI_perm) <- paste0(NI_perm$chr, "_", NI_perm$pos)
        rownames(FLU_perm) <- paste0(FLU_perm$chr, "_", FLU_perm$pos)

        ## for each gene, store each output separately
        perm_NI=NI_perm[rownames(NI_perm) %in% gene , ]
        perm_FLU=FLU_perm[rownames(FLU_perm) %in% gene , ]

        if(j == 1){
          NI_outs <- perm_NI
          FLU_outs <- perm_FLU
        }else{
          NI_outs <- rbind(NI_outs, perm_NI)
          FLU_outs <- rbind(FLU_outs, perm_FLU)
        }
      }, error=function(e){
        #err_genes = c(err_genes, gene)
        #cat("ERROR :",conditionMessage(e), "\n")
      })
    
    ## extract the admixture coefficients
      permuted_adm_NI = NI_outs
      permuted_adm_FLU = FLU_outs
    }

    if(iter==1)
    {
      shuffled_pvals_adm_FLU <-data.frame(x=permuted_adm_FLU$pvals)
      shuffled_pvals_adm_NI <-data.frame(x=permuted_adm_NI$pvals)

      rownames(shuffled_pvals_adm_FLU)=rownames(permuted_adm_FLU)
      rownames(shuffled_pvals_adm_NI)=rownames(permuted_adm_NI)

    } else {
      shuffled_pvals_adm_FLU <- cbind(shuffled_pvals_adm_FLU,x=permuted_adm_FLU$pvals)
      shuffled_pvals_adm_NI <- cbind(shuffled_pvals_adm_NI,x=permuted_adm_NI$pvals)

    }
  }
  system(paste0("mkdir -p ", out_dir, "/permuted_pvalues/"))
  #Writing permutation tests p-values tables
  write.table(shuffled_pvals_adm_FLU,paste0(out_dir, "permuted_pvalues/popDE_FLU.txt"))
  write.table(shuffled_pvals_adm_NI,paste0(out_dir, "permuted_pvalues/popDE_NI.txt"))
}



##########################################################################
### Write results and add additional FDR correction methods to results ###
##########################################################################

## Remove NAs (will be same between real data and permutations)
FLU_results <- FLU_results[complete.cases(FLU_results), ]
NI_results <- NI_results[complete.cases(NI_results), ]

shuffled_pvals_adm_FLU <- shuffled_pvals_adm_FLU[complete.cases(shuffled_pvals_adm_FLU), ]
shuffled_pvals_adm_NI <- shuffled_pvals_adm_NI[complete.cases(shuffled_pvals_adm_NI), ]

dim(FLU_results)[1] == dim(shuffled_pvals_adm_FLU)[1]
dim(NI_results)[1] == dim(shuffled_pvals_adm_NI)[1]

##calculate BH adjusted p-values.
BH_Flu_pvals <- p.adjust(FLU_results$pvals, method= "BH", n=length(FLU_results$pvals))
BH_NI_pvals <- p.adjust(NI_results$pvals, method= "BH", n=length(NI_results$pvals))

##Bind BH adjusted p-values to results table.
results_FLU_BH=cbind(FLU_results, BH_Flu_pvals)
results_NI_BH=cbind(NI_results, BH_NI_pvals)

## Calculate qvalue
q_Flu_pvals <- empPvals(stat=-log10(FLU_results$pvals), stat0=-log10(as.matrix(shuffled_pvals_adm_FLU)), pool = T)
qvals_Flu <- qvalue(p=q_Flu_pvals)$qvalue
q_NI_pvals <- empPvals(stat=-log10(NI_results$pvals), stat0=-log10(as.matrix(shuffled_pvals_adm_NI)), pool = T)
qvals_NI <- qvalue(p=q_NI_pvals)$qvalue

## Bind qvalue to same tables.
results_FLU_BH_q=cbind(results_FLU_BH, qvals_Flu)
results_NI_BH_q=cbind(results_NI_BH, qvals_NI)

## Write results tables
write.table(results_FLU_BH_q, paste0(out_dir, "popDE_FLU.qval.txt"))
write.table(results_NI_BH_q, paste0(out_dir, "popDE_NI.qval.txt"))

############################
## Print sig. PopDE genes ##
############################

#write this as an output and place in summary directory.
sink(paste0(out_dir, "summary_info.txt"))

#what are the PopDE genes in each condition using the different FDR thresholds?
print("Flu PopDE genes with BH <= 0.10")
print(summary(results_FLU_BH$BH_Flu_pvals <= 0.10))
print("NI PopDE genes with BH <= 0.10")
print(summary(results_NI_BH$BH_NI_pvals <= 0.10))
print("Flu PopDE genes with q <= 0.10")
print(summary(results_FLU_BH_q$qvals_Flu <= 0.10))
print("NI PopDE genes with q <= 0.10")
print(summary(results_NI_BH_q$qvals_NI <= 0.10))

#finish printing to file.
sink()
