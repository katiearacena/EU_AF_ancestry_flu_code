##load libraries
library(glmnet)
library(dplyr)
library(stringi)
library(stringr)
library(readxl)

set.seed(917)  # Set seed for reproducibility

## Load command line arguments
args<-commandArgs(TRUE)
## Datatype
data <- args[1]

if(length(args)==0)
{
  print("WARNING: No arguments supplied.")
}
print(args)

## Set directory structure
folder = "2_ancestry_effects_modeling"
## Set directory 3 steps above script.
setwd('../../../')
system(paste0("mkdir -p Outputs/",folder,"/epi_priming/", data, "/"))
out_dir <- paste0("Outputs/",folder,"/epi_priming/", data, "/")

#create output directories
dir.create(file.path(paste0(out_dir, "predicted")))
dir.create(file.path(paste0(out_dir, "weights")))
pred_out = file.path(paste0(out_dir, "predicted"))
weights_out = file.path(paste0(out_dir, "weights"))

## load metadata
metadata <-read.table(paste0("Inputs/metadata/", data,"_metadata.txt"))
## remove unused samples
metadata <- metadata[which(!metadata$Condition=="Mock"),]
metadata <- metadata[which(!metadata$PopDE_set=="0"), ]
metadata <- metadata[metadata$Condition == "NI", ]

## how many samples are included in this analysis (max is 35 inds.)
print(paste0("total number of samples included in analysis = ", dim(metadata)[1]))
   
## load and format RNAseq metadata
RNAseq_metadata <-read.table(paste0("Inputs/metadata/RNAseq_metadata.txt"))
RNAseq_metadata <- RNAseq_metadata[which(!RNAseq_metadata$Condition=="Mock"),]
RNAseq_metadata <- RNAseq_metadata[RNAseq_metadata$Condition == "NI", ]
RNAseq_metadata <- RNAseq_metadata[which(!RNAseq_metadata$PopDE_set=="0"), ]
  
## load FCs
RNAseq_FCs <- read.table(paste0("Outputs/2_ancestry_effects_modeling/PopDR_analysis/RNAseq/results/FC_for_plotting.txt"))
 
## make sure that the dimensions of the reads and metadata match.
dim(RNAseq_metadata)[1]==dim(RNAseq_FCs)[2]
reads_t <- t(RNAseq_FCs)
reads_t <- as.data.frame(reads_t)

## load popDE NI and popDE Flu IFNA scores
GSVA_results <- read.table("Outputs/2_ancestry_effects_modeling/ancestry_score/popDR_GSVEA_SETS_FOR_PLOTTING.txt")
GSVA_results <- GSVA_results[GSVA_results$data == "RNAseq", ]

#############################
## 1. Get IFNA ind. scores ##
#############################

## extract IFNA
IFNA <- GSVA_results[GSVA_results$set=="IFNA", ]
##remove number from rownames and only get first 4 digits
IFNA$Genotyping_ID <- str_sub(rownames(IFNA),1,4)
IFNA_score <- select(IFNA, Genotyping_ID, ind_means)
colnames(IFNA_score)[2] <- "IFNA_score"
metadata <- left_join(metadata, IFNA_score, by="Genotyping_ID")

#############################
## 2. Get IFNG ind. scores ##
#############################

IFNG <- GSVA_results[GSVA_results$set=="IFNG", ]
##remove number from rownames and only get first 4 digits
IFNG$Genotyping_ID <- str_sub(rownames(IFNG),1,4)
IFNG_score <- select(IFNG, Genotyping_ID, ind_means)
colnames(IFNG_score)[2] <- "IFNG_score"
metadata <- left_join(metadata, IFNG_score, by="Genotyping_ID")

#############################
## 3. Get TNFA ind. scores ##
#############################

TNFA <- GSVA_results[GSVA_results$set=="TFNA", ]
##remove number from rownames and only get first 4 digits
TNFA$Genotyping_ID <- str_sub(rownames(TNFA),1,4)
TNFA_score <- select(TNFA, Genotyping_ID, ind_means)
colnames(TNFA_score)[2] <- "TNFA_score"
metadata <- left_join(metadata, TNFA_score, by="Genotyping_ID")

#############################
## 4. Get IL6 ind. scores  ##
#############################

IL6 <- GSVA_results[GSVA_results$set=="IL6", ]
##remove number from rownames and only get first 4 digits
IL6$Genotyping_ID <- str_sub(rownames(IL6),1,4)
IL6_score <- select(IL6, Genotyping_ID, ind_means)
colnames(IL6_score)[2] <- "IL6_score"
metadata <- left_join(metadata, IL6_score, by="Genotyping_ID")

## write metadata with new columns for correlation with prediction later on.
write.table(metadata, paste0(out_dir, "popDR_score_ENR_metadata.txt"))

#######################################################
## load reads (batch and age corrected) for datatype ##
#######################################################

epi_reads <- read.table(paste0("Outputs/2_ancestry_effects_modeling/batch_corrected_cts_matrices/", data, "_batch.age.corrected.txt"))
epi_NI_reads <- epi_reads[, grepl("_NI", colnames(epi_reads))]
colnames(epi_NI_reads) <- str_remove(colnames(epi_NI_reads), "_NI")
cts <- epi_NI_reads

## set alpha parameters to try
alpha=c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
NF <- dim(cts)[2]-1 # (number of samples) - 1
print(paste0("number of samples-1 = ", NF))
models <- c("popDR_RNAseq_IFNA_score", "popDR_RNAseq_IFNG_score", "popDR_RNAseq_TNFA_score", "popDR_RNAseq_IL6_score")

################################
## Elastic net model building ##
################################

for (m in 1:4){
  print(models[m])

  for(i in 1:length(alpha)){
    
    print(alpha[i]) 
    
    #create table to store results
    pred_samp <- data.frame(samp=character(0),ind=character(0),predicted=numeric(0), weight=numeric(0), stringsAsFactors = FALSE)
    
    for(SAMP in 1:dim(cts)[2]){
      
      #################
      ## qqnorm data ##
      #################
      
      ## make sure that metadata and cts match (and remove any samples there is not count information for)
      metadata_i <- metadata[metadata$Genotyping_ID %in% colnames(cts), ]
      colnames(cts)==metadata_i$Genotyping_ID ## make sure == TRUE
      
      #Using alpha
      #Quantile normalize gene exp by column (sample)
      norm_counts<-as.matrix(apply(cts,2,function(x){return(qqnorm(x,plot=F)$x)}))
      
      #Remove test subject(s)
      #SAMP indexes from 1 to n samples
      norm_train<-norm_counts[,-SAMP]
      norm_test<-norm_counts[,SAMP]
      
      #Quantile normalize training samples by row (features)
      trainreads_norm<-as.matrix(apply(norm_train,1,function(x){return(qqnorm(x,plot=F)$x)}))
      
      if (m==1){
        train_model <- metadata$IFNA_score[-SAMP]
      }else if (m==2){
        train_model <- metadata$IFNG_score[-SAMP]
      }else if (m==3){
        train_model <- metadata$TNFA_score[-SAMP]
      }else if (m==4){
        train_model <- metadata$IL6_score[-SAMP]
      }

      #QQ normalize each row in the test sample
      #Note that this could be much more efficient using parallelR (e.g. parSapply)
      
      #Store a new vector for normalizing the test sample
      testreads_normalized<-norm_test
      
      #For each feature
      for (d in 1:length(testreads_normalized)){
        #Define the ECDF (depending on training sample size, this can be replaced with simply the training data -- norm_train)
        a<-ecdf(norm_counts[d,])
        #From this ECDF, return the probability of values being less than the training sample
        probs<-a(norm_test[d])
        #To avoid extreme values outside the range of the training samples, give 0's and 1's a manageable quantile (e.g. 0.99 and .01)
        #Note depending on the sample size, consider changing this number (e.g. to 1/N and 1-1/N respectively)
        probs[probs==1]<-.99
        probs[probs==0]<-.01
        #Given this probability, return the associated value from a standard normal that falls into the same quantile
        testreads_normalized[d]<-qnorm(probs)
      }
      
      ################################
      ## Elastic-net model building ##
      ################################
      
      #Using N-fold internal CV, train the elastic net model using the training data
      #Note with larger sample sizes, N-fold internal CV becomes intractable
      model<-cv.glmnet(trainreads_norm,train_model,nfolds=NF,alpha=alpha[i],standardize=F)
      
      #Predict Ab using the test sample from parameters that minimized MSE during internal CV
      predicted_SAMP<-predict(model,newx=t(testreads_normalized),s="lambda.min")
      
      #Extract weights for this model
      weights_SAMP<-unlist(coef(model,lambda="lambda.min"))[,1]
      
      #Write out results for later concatenation
      write.table(predicted_SAMP,paste0(out_dir, "predicted/",models[m], "_predicted_QQ_n",NF,"_a",alpha[i],"_s",SAMP,".txt"),quote=F,row.names=F,col.names=F)
      write.table(weights_SAMP,paste0(out_dir, "weights/", models[m], "_weights_QQ_n",NF,"_a",alpha[i],"_s",SAMP,".txt"),quote=F,row.names=F,col.names=F)
      
      pred_samp[nrow(pred_samp)+1,] <- c(SAMP,(as.character(metadata[SAMP,"Genotyping_ID"])),predicted_SAMP, weights_SAMP)
    }
    
    write.table(pred_samp, paste0(out_dir, models[m], "_predicted_",data,"_n",NF,"_a",alpha[i],".txt"), quote = FALSE, sep = ",", row.names = F)
  }
}
