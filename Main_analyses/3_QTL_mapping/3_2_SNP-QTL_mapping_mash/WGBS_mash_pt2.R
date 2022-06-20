## part 1 of the WGBS mash pipeline extracts and saves the data in the proper format to run mash.
## part 2 of the pipeline applies mash to the data and saves results.

# Load libraries
library(ashr)
library(mashr)
library(gplots)
library(viridis)
library(ggplot2)
library(corrplot)
library(rmeta)
library(dplyr)
set.seed(321)

############################
## Load command line args ##
############################

args<-commandArgs(TRUE)

## Declare data type
part <- args[1]


if(length(args)==0)
  {
    print("WARNING: No arguments supplied.")
  }
print(args)

part <- as.character(part)
data <- "WGBS"

##############################################
## Create directory structure and load data ##
##############################################

folder = "3_QTL_mapping"
## Set directory 3 steps above script.
setwd('../../../')

## Create directory structure to save outputs.
system(paste0("mkdir -p Outputs/",folder,"/SNP-QTL_mash/", data, "/part", part, "_outputs/"))
out_dir <- paste0("Outputs/",folder,"/SNP-QTL_mash/", data, "/part", part, "_outputs/")

load(paste0("Outputs/",folder,"/SNP-QTL_mash/WGBS/mash_random.Rdata"))
x <- load(paste0("Outputs/",folder,"/SNP-QTL_mash/WGBS/data_STRONG_p", part, ".Rdata"))
data_STRONG = get(x)

# Remove the old object since you've stored it in x
rm(x)

m2 = mash(data_STRONG, g = get_fitted_g(m), fixg = TRUE)

#################################
## Extract Posterior Summaries ##
#################################

# Get effects that are “significant”, which  means they have lfsr less than .10 in at least one condition
thresh <- .10
m.pairwise_PM <- get_pairwise_sharing(m2, lfsr_thresh=thresh, factor = 0.5)
m.sig <- get_significant_results(m2, thresh=thresh, sig_fn=get_lfdr)

################################################
## Identify condition-specific and shared QTL ##
################################################

#write this as an output and place in summary directory. ***
sink(paste0(out_dir, "summary_info.txt"))

print("how many are significant in at least one condition")
print(length(m.sig))
print("how many are significant in just Flu")
length(get_significant_results(m2, thresh=thresh, sig_fn=get_lfsr, conditions=1))
print("how many are significant in just NI")
length(get_significant_results(m2, thresh=thresh, sig_fn=get_lfsr, conditions=2))

#finish printing to file.
sink()

## write posterior outs
write.table(get_lfsr(m2), paste0(out_dir, "lfsr_output.txt"), quote = FALSE)
write.table(get_lfdr(m2), paste0(out_dir, "lfdr_output.txt"), quote = FALSE)
write.table(get_pm(m2), paste0(out_dir, "posteriorMeans.txt"), quote = FALSE)
write.table(get_psd(m2), paste0(out_dir, "posteriorStandardDevs.txt"), quote = FALSE)

## save R object
save(m2, file=paste0(out_dir, "mash_output.Rdata"))