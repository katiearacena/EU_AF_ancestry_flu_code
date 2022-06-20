### script to combine all the null chromosomes (the outputs of QTL_integ_footprinting_NULL)

## load inputs 
folder = "3_QTL_mapping"
## Set directory 3 steps above script.
setwd('../../../')

out_dir <- paste0("Outputs/",folder,"/QTL_integ_TF_enrichments/")

## for NI
NI_dir <- paste0(out_dir, "NI")
NI_names <- list.files(NI_dir, "NULL.", all.files=TRUE, full.names=TRUE)
NI_files <- lapply(NI_names, readRDS)
##make sure all 22 chromosomes are present
length(NI_files)==22

NI.NULL.ALL.RDS <- do.call(rbind.data.frame, NI_files)
dim(NI.NULL.ALL.RDS)
saveRDS(NI.NULL.ALL.RDS, paste0(NI_dir, "/NULL.ALL.RDS"))

## for Flu
Flu_dir <- paste0(out_dir, "Flu")
Flu_names <- list.files(Flu_dir, "NULL.", all.files=TRUE, full.names=TRUE)
Flu_files <- lapply(Flu_names, readRDS)
##make sure all 22 chromosomes are present
length(Flu_files)==22

Flu.NULL.ALL.RDS <- do.call(rbind.data.frame, Flu_files)
dim(Flu.NULL.ALL.RDS)
saveRDS(Flu.NULL.ALL.RDS, paste0(Flu_dir, "/NULL.ALL.RDS"))
