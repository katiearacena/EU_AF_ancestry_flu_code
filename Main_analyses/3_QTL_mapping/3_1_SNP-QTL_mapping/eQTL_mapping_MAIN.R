## Load required libraries.
library(MatrixEQTL)
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
## Declare directory for temp files
temp_dir <- args[4]
## Declare output directory.
res_dir <- args[5]

if(length(args)==0)
  {
    print("WARNING: No arguments supplied.")
  }
print(args)


# Perform 10 permutations to use for FDR correction.
iterations <- 10
print(data)


#############################################################
## CREATE DIRECTORY STRUCTURE, SET FUNCTIONS, & LOAD DATA  ##
#############################################################

folder = "3_QTL_mapping"
## Set directory 3 steps above script.
setwd('../../../')

## Create directory structure to save outputs.
system(paste0("mkdir -p Outputs/",folder,"/SNP-QTL_mapping/", data, "/", condition, "/", res_dir, "/"))
out_dir <- paste0("Outputs/",folder,"/SNP-QTL_mapping/", data, "/", condition, "/", res_dir, "/")

## Create temp directory.
system(paste0("mkdir -p ", out_dir, "temp_files/"))

## Create directory for raw results.
system(paste0("mkdir -p ", out_dir, "/raw_results/"))

## Create directory for qqplots.
system(paste0("mkdir -p ", out_dir, "/qqplots/"))

## Create directory for pvalue_histograms (if plan to create them).
system(paste0("mkdir -p ", out_dir, "/pval_hists"))

## Erase any previous files from and re-create the temp file directory.
system(paste0("rm -rf ",out_dir, "temp_files/", temp_dir, "/"))
system(paste0("mkdir -p ",out_dir, "temp_files/", temp_dir, "/"))

##Load metadata
metadata_whole = read.table(paste0("Inputs/metadata/", data, "_metadata.txt"))
## Load age and batch corrected voom normalized read counts
reads_whole = read.table(paste0("Outputs/2_ancestry_effects_modeling/batch_corrected_cts_matrices/",data, "_batch.age.corrected.txt"), header=TRUE, row.names=1)

## Load positions file
genepos = read.table(paste0("Inputs/QTL_mapping/",data, "_positions.txt"),header = TRUE, stringsAsFactors = FALSE)
## Load SNP positions
snpspos = read.table(paste0("Inputs/QTL_mapping/SNP_positions.txt"),header = TRUE, stringsAsFactors = FALSE)
## Load genotype info for SNPs
gtypes = read.table(paste0("Inputs/QTL_mapping/SNP_genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)

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
## Build matrixEQTL input: expression tables: regressing out first n PCs ##
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
    expression = pca_rm(reads, pc_set)
}

## Now have expression matrix with normalized reads and n PCs regressed. 


#########################################################################################
## Build matrixEQTL input: covariates tables: transpose metadata & add 1st genotype PC ##
#########################################################################################

## Perform genotype PC analysis and clean genotype data.
## Remove condition from rownames of metadata_individuals.
metadata_individuals=metadata_whole[which(!duplicated(metadata_whole$Genotyping_ID)),]
metadata_individuals=metadata_individuals[order(metadata_individuals$Genotyping_ID),]
rownames(metadata_individuals)=metadata_individuals$Genotyping_ID

## Set the column names of gtypes to be the samples.
samples=colnames(gtypes)
gtypes_pca=data.frame(snp_id=rownames(gtypes),gtypes)

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

## Check to make sure that rownames are correct and match gene_pos.
expression=expression[which(rownames(expression) %in% genepos$Gene_ID),]

##########################################################################################
## Before calling matrixEQTL subset individuals present in the condition to analyze and ##
## remove genes and SNPs from genotype and expression tables                            ##
## for which there is no available position (These wouldn't be tested for cis-EQTL,     ##
## but had useful info to include in PC analyses performed above)                       ##
##########################################################################################

## Subset the genotypes file with the individuals present in each condition.
genotypes=gtypes[,which(colnames(gtypes) %in% colnames(covariates))]
genotypes=genotypes[which(rownames(genotypes) %in% snpspos$snp),]
genotypes=genotypes[,order(colnames(genotypes))]
length(which(rownames(genotypes)!=snpspos$snp))

## Check input data files congruence 

## Check for sample congruence. These should all be 0. If it is not 0 there is an issue. 
length(which(rownames(metadata)!=colnames(covariates)))
length(which(rownames(metadata)!=colnames(expression)))
length(which(rownames(metadata)!=colnames(genotypes)))
length(which(rownames(genotypes)!=snpspos$snp))

## Gene-wise check. This step removes any genes for which there is not expression data. 
#use negate now since genepos and expression are not in same order.
`%nin%` = Negate(`%in%`)

genepos_trimmed <- genepos[(genepos$Gene_ID %in% rownames(expression)),]
genepos <- genepos_trimmed
length(which(rownames(expression) %nin% genepos_trimmed$Gene_ID))
length(which(rownames(expression) %nin% genepos$Gene_ID))
## If all 0s, everything is coherent.

length(genepos_trimmed$Gene_ID)
length(genepos$Gene_ID)
length(rownames(expression))
## This output will tell you how many genes there are expression data for. They should all have the same length.

#######################################
## Save matrixEQTL temp input files ###
#######################################

## Save files.
snps_positions_file_name="Inputs/QTL_mapping/SNP_positions.txt"
gene_positions_file_name=paste0("Inputs/QTL_mapping/", data, "_positions.txt")
expression_file_name=paste0(out_dir, "temp_files/", temp_dir,"/expression.txt")
covariates_file_name=paste0(out_dir, "temp_files/", temp_dir,"/covariates.txt")
SNP_file_name=paste0(out_dir, "temp_files/", temp_dir,"/genotypes.txt")

## Write files.
write.table(genotypes,SNP_file_name, quote=F, sep="\t", row.names=TRUE)
write.table(expression,expression_file_name, quote=F, sep="\t", row.names=TRUE)
write.table(covariates,covariates_file_name, quote=F, sep="\t", row.names=TRUE)

## In this loop iter=0 runs the actual analyses, iters 1 to iterations (10), run permutations for FDR correction.

permuted_pvalues_folder=paste0(out_dir,"/raw_results/")
for(iteration in 0:iterations){
    
  ###################################################
  ## Permute genotype data (only for iterations>0) ##
  ###################################################
  
  if(iteration>0){
    
    cols<-colnames(genotypes)
    cols.perm<-sample(cols)
    if(iteration==1){
      random_individuals_df=data.frame(cols.perm)
    }else{
      random_individuals_df=cbind(random_individuals_df,cols.perm)
    }
    genotypes<-genotypes[,cols.perm]
    colnames(genotypes)<-cols
    write.table(genotypes,SNP_file_name, sep="\t", quote = FALSE)
  }

    ###############################
    ## Prepare & Run Matrix EQTL ##
    ###############################
    
    ## Load phenotype data.
    gene = SlicedData$new();
    gene$fileDelimiter = "\t";     # the TAB character
    gene$fileOmitCharacters = "NA"; # denote missing values;
    gene$fileSkipRows = 1;          # one row of column labels
    gene$fileSkipColumns = 1;       # one column of row labels
    gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    gene$LoadFile(expression_file_name);
    
    ## Load covariates
    ## NOTE: If you are not regressing out any covariates add this line.
    # covariates_file_name = character()

    cvrt = SlicedData$new();
    cvrt$fileDelimiter = "\t";      # the TAB character
    cvrt$fileOmitCharacters = "NA"; # denote missing values;
    cvrt$fileSkipRows = 1;          # one row of column labels
    cvrt$fileSkipColumns = 1;       # one column of row labels
    if(length(covariates_file_name)>0) {
        cvrt$LoadFile(covariates_file_name);
    }
    
    ## Load genotype data
    
    snps = SlicedData$new();
    snps$fileDelimiter = "\t";      # the TAB character
    snps$fileOmitCharacters = "NA";
    snps$fileOmitCharacters = "-9" # denote missing values;
    snps$fileSkipRows = 1;          # one row of column labels
    snps$fileSkipColumns = 1;       # one column of row labels
    snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    snps$LoadFile(SNP_file_name)
    
    ## Set MatrixeQTL options

    useModel = modelLINEAR
    output_file_name_cis = tempfile()
    pvOutputThreshold_cis = 1
    pvOutputThreshold = 0;
    errorCovariance = numeric()
    cisDist = 1e5
    output_file_name = tempfile()
    output_file_name_cis = tempfile()
    
    ## Begin MatrixeQTL

    me = Matrix_eQTL_main(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = output_file_name,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = TRUE,
        output_file_name.cis = output_file_name_cis,
        pvOutputThreshold = pvOutputThreshold,
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos = snpspos,
        genepos = genepos,
        cisDist = cisDist,
        pvalue.hist = "qqplot",
        min.pv.by.genesnp = TRUE,
        noFDRsaveMemory = FALSE);

    ## Create and save qqplot for each iteration.
    png(filename = paste0(out_dir, "qqplots/qqplot", iteration,".png"))
    plot(me, pch = 16, cex = 0.7)
    dev.off()

 ## Can also create pvalue histograms (but takes extra time).
  
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold = pvOutputThreshold,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = 100,
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = FALSE);
  
 ## Create pvalue histogram
   png(filename = paste0(out_dir, "/pval_hists/pval_histogram" ,iteration,".png"))
   plot(me, col="grey")
   dev.off()
    
##############################################################
#### 9. Write temporal output files ##########################
##############################################################
    
    unlink(output_file_name_cis);

    if(iteration==0){
        write.table(me$cis$eqtls, file = paste0(out_dir,"/raw_results/result_original.txt"))}else{
            write.table(me$cis$eqtls, file = paste0(out_dir,"/raw_results/result_permuted_",iteration,".txt"))}
}


#######################################################################################################
## Write Best SNP-gene associations files for Delta-PVE analyses & prepare files for FDR corrections ##
#######################################################################################################

## Select the top cis-SNP for each gene in true and permuted files.
for(iteration in 0:iterations)
{ 
    if(iteration==0){
        event=read.table(paste0(out_dir,"/raw_results/result_original.txt"),header=TRUE)
        event.sort<-event[order(event[,2],event[,4]),]
        event.bestQTL<-event.sort[!duplicated(event.sort$gene),]
        event.bestQTL<-event.bestQTL[order(event.bestQTL[,4]),]
    }else{
        event=read.table(paste0(permuted_pvalues_folder,"/result_permuted_",iteration,".txt"),header=TRUE)
        event.sort<-event[order(event[,2],event[,4]),]
        event.bestQTL<-event.sort[!duplicated(event.sort$gene),]
        event.bestQTL<-event.bestQTL[order(event.bestQTL[,5]),]
    }
    if(iteration==0){
        original_best_EQTL=event.bestQTL
    }else{
        if(iteration==1)
        {   
            permuted1_best_EQTL=event.bestQTL
            #only pull out the p-value.
            Permutation_Input = event.bestQTL[4]
        }else{
            Permutation_Input=cbind(Permutation_Input,event.bestQTL[4])}
    }
}


#########################################
## Use qvalue package to calculate FDR ##
#########################################

## Calculate qvalues using the best_SNPs files.
emp_pvalues <- empPvals(stat=-log10(original_best_EQTL[,4]), stat0=-log10(as.matrix(Permutation_Input)), pool = T)
qvalues <- qvalue(p=emp_pvalues)$qvalue

## Bind qvalue output with results_best_SNPs file (assumes order is the same in both).
results_best_SNPs_with_qval <- cbind(original_best_EQTL, qvalues)
## Write table with this column added.
write.table(x=results_best_SNPs_with_qval, file= paste0(out_dir,"results_best_SNPs_with_qval.txt"), quote=FALSE)

## Save permuted best SNPs for 1 permutation
write.table(permuted1_best_EQTL, file = paste0(out_dir,"permuted1_best_SNPs.txt"), quote=FALSE)

####################################################
## Create qqplot using p-values of best SNPs only ##
####################################################

## Make and save best SNPs qqplot.
png(filename =  paste0(out_dir,"best_SNPs_qqplot.png"))
qqplot(-log10(permuted1_best_EQTL[,4]), -log10(results_best_SNPs_with_qval[,4]), ylim=  c(0, 40), xlim= c(0, 40), main ="Best SNPs qqplot", xlab = "-log10(theoretical p-values)", ylab = "-log10(observed p-values)")
abline(c(0,1),col="red")
dev.off()

########################################
## Add SNP ID to results_original.txt ##
########################################

rsIDs <- read.table("Inputs/ref/rsIDs.for.matrixeqtl.txt")
rsID_snps_only <- dplyr::select(rsIDs, rsID, snps)
result_original <- read.table(paste0(out_dir,"/raw_results/result_original.txt"),header=TRUE)
result_original_with_rsID <- right_join(rsID_snps_only, result_original, by = "snps")
write.table(result_original_with_rsID, paste0(out_dir,"/raw_results/result_original_rsIDs.txt"))