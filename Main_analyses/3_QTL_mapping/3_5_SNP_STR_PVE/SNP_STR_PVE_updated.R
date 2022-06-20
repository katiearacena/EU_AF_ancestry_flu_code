#######################################################################################

.libPaths("/project2/lbarreiro/users/Onta/R_3.6.3-no-openblas/")
## Load required libraries.
library(edgeR)
library(limma)
library(gdsfmt)
library(SNPRelate)
library(qvalue)
library(gplots)
library(viridis)
library(ggplot2)
library(corrplot)
library(rmeta)
library(dplyr)
library(bsseq)
library(data.table)
library(relaimpo)

#################################
##### LOAD COMMAND LINE ARGS ####
#################################

args<-commandArgs(TRUE)
## Declare condition equal to "NI or "FLU" to obtain cis-EQTLs for each condition.
condition <- args[1]
## Declare data type.
data <- args[2]
## Declare number of expression PCs to regress out.
expPCs_reg <- args[3]


if(length(args)==0)
{
  print("WARNING: No arguments supplied.")
}
print(args)


print(data)

CONDITION <- toupper(condition)
#############################################################
## CREATE DIRECTORY STRUCTURE, SET FUNCTIONS, & LOAD DATA  ##
#############################################################

folder = "3_QTL_mapping"
## Set directory 3 steps above script.
setwd('../../../')


## Load genotype info for STRs and SNPs
STRgtypes = read.table(paste0("Inputs/QTL_mapping/STR_genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)
SNPgtypes = read.table(paste0("Inputs/QTL_mapping/SNP_genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)

## Create directory structure to save outputs.
system(paste0("mkdir -p Outputs/",folder,"/PVE_expression_matched_MEQTL/", data,"/"))
out_dir <- paste0("Outputs/",folder,"/PVE_expression_matched_MEQTL/", data,"/")

##Load metadata
metadata_whole = read.table(paste0("Inputs/metadata/", data, "_metadata.txt"))

## load counts

if (data=="WGBS"){
    load("Inputs/counts_matrices/WGBS_filtered.counts.BSobj.RData")
    BS.cov <- getCoverage(BSobj.fit.nolow)

    ## Filter based on coverage
    #To minimize noise in methylation estimates due to low-coverage data, we restricted analysis to CpG sites
    # with coverage of ≥4 sequence reads in at least half of the samples in each condition. DSS accounts for low coverage
    # in the model so it is not necessary before.

    keepLoci.ex <- which(rowSums(BS.cov[, BSobj.fit.nolow$Condition == "Flu"] >= 4) > 17 &
                       rowSums(BS.cov[, BSobj.fit.nolow$Condition == "NI"] >= 4) > 17 )
    length(keepLoci.ex)
    BSobj.fit.nolow <- BSobj.fit.nolow[keepLoci.ex,]
    
    ## get cts matrix from bsseq object
    reads_whole <- getMeth(BSobj.fit.nolow, type="raw")
    
    reads_whole=reads_whole[,which(colnames(reads_whole) %in% rownames(metadata_whole))]
    # Make sure that the samples are ordered.
    reads_whole=reads_whole[order(rownames(reads_whole)),order(colnames(reads_whole))]

    ## This is a final check to make sure the order of elements in the metadata_whole and reads_whole.
    ## If everything is correct the length = 0.
    length(which(rownames(metadata_whole)!=colnames(reads_whole)))
    
    ## remove sites with no variation
    no.var <- rownames(reads_whole[rowVars(reads_whole, na.rm=TRUE) == 0, ])
    reads_whole.lim <- subset(reads_whole, !(rownames(reads_whole) %in% no.var))

    ## quantile normalize
    quantile_expression <- matrix(, nrow = nrow(reads_whole.lim), ncol = ncol(reads_whole.lim))

    for (j in 1:nrow(reads_whole.lim )){
        exp <- reads_whole.lim [j,]
        exp_QN <- qqnorm(exp, plot.it =F)$x
        quantile_expression[j,] <- exp_QN
    }

    colnames(quantile_expression) <- colnames(reads_whole.lim)
    rownames(quantile_expression) <- rownames(reads_whole.lim)
    reads_whole <- quantile_expression
    
    } else {
    
    ## Load age and batch corrected voom normalized read counts
    reads_whole = read.table(paste0("Outputs/2_ancestry_effects_modeling/batch_corrected_cts_matrices/",data, "_batch.age.corrected.txt"), header=TRUE, row.names=1)
}


## Remove string from chipseq rownmames to make sure they match counts matrix
rownames(metadata_whole)<- sub("_H3K27ac", "", rownames(metadata_whole))
rownames(metadata_whole)<- sub("_H3K4me1", "", rownames(metadata_whole))
rownames(metadata_whole)<- sub("_H3K27me3", "", rownames(metadata_whole))
rownames(metadata_whole)<- sub("_H3K4me3", "", rownames(metadata_whole))

# EQTL_set exists as an extra check on which samples to include in EQTL analysis.
metadata_whole=metadata_whole[which(metadata_whole$EQTL_set==1),]
## Select only samples for which genotype data is available in the global metadata table, and order samples:
metadata_whole=metadata_whole[order(rownames(metadata_whole)),]
## Subset the corresponding columns in the reads matrix, and order samples and genes.
#this removes samples that are in the reads_whole matrix that are not in the metadata_whole matrix.
reads_whole=reads_whole[,which(colnames(reads_whole) %in% rownames(metadata_whole))]
# Make sure that the samples are ordered.
reads_whole=reads_whole[order(rownames(reads_whole)),order(colnames(reads_whole))]

## This is a final check to make sure the order of elements in the metadata_whole and reads_whole.
## If everything is correct the length = 0.
length(which(rownames(metadata_whole)!=colnames(reads_whole)))


##################################################
## Filter expression and metadata per condition ##
##################################################


## Filter metadata based on condition.
metadata=metadata_whole[which(metadata_whole$Condition==condition),]

## Clean factor variables, and mean center numeric ones.
metadata$Condition=factor(metadata$Condition)
metadata$Genotyping_ID=factor(metadata$Genotyping_ID)
metadata$Batch=factor(metadata$Batch)

## Filter expression per condition
reads=reads_whole[,which(colnames(reads_whole) %in% rownames(metadata))]
## Check again the coherence of samples order after filtering for condition. Should be 0.
length(which(colnames(reads)!=rownames(metadata)))

## Shift from sampleIDs to Genotyping_IDs (there only will be one sample per genotype in every analysis).
## This removes the condition from the col names of  reads and row names of metadata.
colnames(reads)=metadata$Genotyping_ID
rownames(metadata)=metadata$Genotyping_ID

## Recover alphabetical order with the new IDs.
reads=reads[,order(colnames(reads))]
metadata=metadata[order(metadata$Genotyping_ID),]
## Check to make sure length = 0.
length(which(colnames( reads)!=rownames(metadata)))


###########################################################################
## Build model input: expression tables: regressing out first n PCs ##
###########################################################################

if(expPCs_reg==0){
  expression=reads
}else{
  ## Empirically remove PCs from the phenotype data.
  pc_set=c(1:expPCs_reg)

  ## Regress those out.
  pca_rm <- function(input_data, pc_set) {
    pca = prcomp(t(input_data), na.action = na.omit)
    new = input_data
    new = apply(new, 1, FUN = function(x){return(lm(as.numeric(x) ~ -1 + pca$x[, as.numeric(pc_set)])$resid)})
    new = t(new)
    colnames(new) = colnames(input_data)
    rownames(new) = rownames(input_data)
    return(new)
  }
  expression = pca_rm(na.omit(reads), pc_set)
}

if(data == "WGBS"){
    saveRDS(expression, paste0(out_dir, condition, "_expPC_expression.RDS"))
    }
## Now have expression matrix with normalized reads and n PCs regressed.

## expression <- readRDS(paste0(out_dir, condition, "_expPC_expression.RDS"))

#########################################################################################
## Build input: covariates tables: transpose metadata & add 1st genotype PC ##
#########################################################################################

## Perform genotype PC analysis and clean genotype data.
## Remove condition from rownames of metadata_individuals.
metadata_individuals=metadata_whole[which(!duplicated(metadata_whole$Genotyping_ID)),]
metadata_individuals=metadata_individuals[order(metadata_individuals$Genotyping_ID),]
rownames(metadata_individuals)=metadata_individuals$Genotyping_ID

## If the covariate table with genotype PC1 was not generated for SNP analysis,
## perform PC analysis using the SNP genotypes
if(file.exists(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/",condition,"/",data,"_",condition,"/temp_files/",data,"_",condition,"/covariates.txt"))){
  covariates=read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/",condition,"/",data,"_",condition,"/temp_files/",data,"_",condition,"/covariates.txt"))
  names(covariates) <- gsub("AF22", "EU22", names(covariates))
  names(covariates) <- gsub("AF36", "EU36", names(covariates))
  names(covariates) <- gsub("AF38", "EU38", names(covariates))
  covariates <- covariates[order(names(covariates))]
}else{
  ## Set the column names of gtypes to be the samples.
  SNPgtypes=read.table(paste0("Inputs/QTL_mapping/SNP_genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)
  samples=colnames(SNPgtypes)
  gtypes_pca=data.frame(snp_id=rownames(SNPgtypes),SNPgtypes)

  ## NOTE: IF YOU RUN INTO AN ERROR: RUN snpgdsClose(genofile) TO RESET FILE OPEN/CLOSE
  # snpgdsClose(genofile)

  ## This creates a SNP genotype dataset from the gtypes_pca matrix.
  snpgdsCreateGeno(paste0(out_dir, "temp_files/", temp_dir, "/GDS_genotypes.gds"),
                   genmat = as.matrix(gtypes_pca[, samples]),
                   sample.id = unique(samples),
                   snp.id = gtypes_pca$snp_id,
                   snpfirstdim=TRUE)

  ## This command tells you the total number of samples and SNPs in the .gds file.
  snpgdsSummary(paste0(out_dir, "temp_files/", temp_dir, "/GDS_genotypes.gds"))

  ## Load .gds file in as genofile.
  genofile <- snpgdsOpen(paste0(out_dir, "temp_files/", temp_dir, "/GDS_genotypes.gds"))

  ## Perform a PCA on genofile/ genotype information.
  pca <- snpgdsPCA(genofile)
  ## Subset the first PC
  tab <- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[,1],    # the first eigenvector
                    stringsAsFactors = FALSE)

  ## Create covariates table.
  ## Make sure that the pcs_genotypes file is properly labeled and ordered.
  pcs_genotypes=tab[which(tab$sample.id %in% rownames(metadata)),]
  pcs_genotypes=pcs_genotypes[order(pcs_genotypes$sample.id),]
  length(which(rownames(metadata) !=pcs_genotypes$sample.id))
  metadata$PC1=pcs_genotypes$PC1

  # batch and age are already accounted for in input reads.
  # PC1 of genotype data is the only covariate included (to account for population structure).
  covariates=t(model.matrix(~PC1,data=metadata))
  covariates=t(matrix((covariates[2:nrow(covariates),])))
  colnames(covariates) <- rownames(metadata)
}



##########################################################################################
## Get the best SNPs and STRs                                                           ##
## filter the genes and variant from genotype and expression tables                     ##
## to reduce the test number                                                            ##
##########################################################################################
library(dplyr)
## get lists of the union/intersect of genes popDE SNPs and popDE STRs.
bestSNPs <- read.table(paste0("Outputs/3_QTL_mapping/SNP-QTL_mapping/",data,"/", condition,"/", data,"_", condition,"/results_best_SNPs_with_qval.txt"))
#rownames(bestSNPs) <- bestSNPs$V1; bestSNPs$V1 <- NULL
bestSNPs <- arrange(bestSNPs, gene)

bestSTRs <- read.table(paste0("Outputs/3_QTL_mapping/STR-QTL_mapping/",data,"/", condition,"/", data,"_", condition,"/results_best_STRs_with_qval.txt"))
#rownames(bestSTRs) <- bestSTRs$V1; bestSTRs$V1 <- NULL
bestSTRs <- arrange(bestSTRs, gene)

genes_of_interests <- union(bestSNPs$gene, bestSTRs$gene)

## extract the genotypes of the significant SNPs
## only keep snp and gene columns
#bestSNPs_subset <- select(bestSNPs, snps, gene)
bestSNPs_subset <- bestSNPs[c("snps", "gene")]
## subset genotypes to pertinent ones
SNPgtypes$snp <- row.names(SNPgtypes)
SNPgenotypes_subset <- right_join(SNPgtypes, bestSNPs_subset, by = c("snp"="snps"))
rownames(SNPgenotypes_subset) <- SNPgenotypes_subset$gene; SNPgenotypes_subset$gene <- NULL; SNPgenotypes_subset$snp <- NULL
## Change missing data from -9 to NA.
SNPgenotypes_subset[SNPgenotypes_subset == -9] <- NA

## extract the genotypes of the significant STRs
## only keep str and gene columns
#bestSTRs_subset <- select(bestSTRs, snps, gene)
bestSTRs_subset <- bestSTRs[c("snps", "gene")]
## subset genotypes to pertinent ones
STRgtypes$str <- row.names(STRgtypes)
STRgenotypes_subset <- right_join(STRgtypes, bestSTRs_subset, by = c("str" = "snps"))
rownames(STRgenotypes_subset) <- STRgenotypes_subset$gene; STRgenotypes_subset$gene <- NULL; STRgenotypes_subset$str <- NULL
## Change missing data from -9 to NA.
STRgenotypes_subset[STRgenotypes_subset == -9] <- NA


## Subset the genotypes file with the individuals present in each condition.
SNPgenotypes_subset <- SNPgenotypes_subset[,which(colnames(SNPgenotypes_subset) %in% colnames(covariates))]
SNPgenotypes_subset <- SNPgenotypes_subset[,order(colnames(SNPgenotypes_subset))]

STRgenotypes_subset <- STRgenotypes_subset[,which(colnames(STRgenotypes_subset) %in% colnames(covariates))]
STRgenotypes_subset <- STRgenotypes_subset[,order(colnames(STRgenotypes_subset))]



## Subset and order the expression and genotype file
EXPRESSION <- expression[rownames(expression) %in% genes_of_interests,]

genes_of_interests <- genes_of_interests[genes_of_interests %in% rownames(EXPRESSION)]

SNPGENOTYPES_SUBSET <- SNPgenotypes_subset[genes_of_interests,]
STRGENOTYPES_SUBSET <- STRgenotypes_subset[genes_of_interests,]


#############################################################################################
## Model expression of each gene with genotype PC1, best SNP genotype, best STR genotype   ##
#############################################################################################

output = as.data.frame(matrix(nrow = length(genes_of_interests), ncol = 3))
names(output) = c("PC1", "SNP", "STR")
row.names(output) = genes_of_interests

for(j in 1:length(genes_of_interests)){
  tryCatch({
    if(j%%10 == 0) print(j)

    gene <- row.names(output)[j]
    COVARIATES <- as.data.frame(t(covariates))
    names(COVARIATES) <- "PC1"
    if (gene %in% row.names(SNPGENOTYPES_SUBSET)){
        COVARIATES$SNP <- as.numeric(SNPGENOTYPES_SUBSET[gene, rownames(COVARIATES)])
    }
    if (gene %in% row.names(STRGENOTYPES_SUBSET)){
        COVARIATES$STR <- as.numeric(STRGENOTYPES_SUBSET[gene, rownames(COVARIATES)])
    }
    vars <- names(COVARIATES)
    
    # relaimpo
    gene_model=lm(reformulate(termlabels = vars, response = 'EXPRESSION[gene,]'), data=COVARIATES)
    rel_impo=suppressWarnings(calc.relimp(gene_model,type="lmg"))
    for(var in vars){
      output[gene,var] = rel_impo$lmg[var]
    }
  }, error=function(e){
    #err_genes = c(err_genes, gene)
    #cat("ERROR :",conditionMessage(e), "\n")
  })  
}
 
write.table(output, paste0(out_dir, data,"_",condition,"_PVE.txt"), quote = FALSE)





