
susie_dir <- "/project2/xinhe/kevinluo/flu_regulatoryQTLs/susie"

#### Distance prior
susie_prior <- "dist"

cat("Combine susie results for WGBS Flu ... \n")
SusieResult_combined <- data.frame()
phenotype_list <-  paste0("WGBS_Flu_full_E", 3:7)

for( phenotype_name in phenotype_list) {

  ## Load susie input data
  cat("Load", phenotype_name, "results ... \n")
  ## checking genes missing susie result
  susie_output_dir <- paste0(susie_dir, "/output/", phenotype_name, "/susie_", susie_prior, "_prior")
  SusieResult <- readRDS(paste0(susie_output_dir, "/susie_", phenotype_name, "_", susie_prior, "_prior_sumstats.rds"))
  SusieResult_combined <- rbind(SusieResult_combined, SusieResult)
}

phenotype_name <- "WGBS_Flu"
susie_output_dir <- paste0(susie_dir, "/output/", phenotype_name, "/susie_", susie_prior, "_prior")
dir.create(susie_output_dir, showWarnings = F, recursive = T)
saveRDS(SusieResult_combined, paste0(susie_output_dir, "/susie_", phenotype_name, "_", susie_prior, "_prior_sumstats.rds"))


cat("Combine susie results for WGBS NI ... \n")
SusieResult_combined <- data.frame()
phenotype_list <-  paste0("WGBS_NI_full_E", 3:7)

for( phenotype_name in phenotype_list) {

  ## Load susie input data
  cat("Load", phenotype_name, "results ... \n")
  ## checking genes missing susie result
  susie_output_dir <- paste0(susie_dir, "/output/", phenotype_name, "/susie_", susie_prior, "_prior")
  SusieResult <- readRDS(paste0(susie_output_dir, "/susie_", phenotype_name, "_", susie_prior, "_prior_sumstats.rds"))
  SusieResult_combined <- rbind(SusieResult_combined, SusieResult)
}

phenotype_name <- "WGBS_NI"
susie_output_dir <- paste0(susie_dir, "/output/", phenotype_name, "/susie_", susie_prior, "_prior")
dir.create(susie_output_dir, showWarnings = F, recursive = T)
saveRDS(SusieResult_combined, paste0(susie_output_dir, "/susie_", phenotype_name, "_", susie_prior, "_prior_sumstats.rds"))

#### Uniform prior
susie_prior <- "uniform"

cat("Combine susie results for WGBS Flu ... \n")
SusieResult_combined <- data.frame()
phenotype_list <-  paste0("WGBS_Flu_full_E", 3:7)

for( phenotype_name in phenotype_list) {

  ## Load susie input data
  cat("Load", phenotype_name, "results ... \n")
  ## checking genes missing susie result
  susie_output_dir <- paste0(susie_dir, "/output/", phenotype_name, "/susie_", susie_prior, "_prior")
  SusieResult <- readRDS(paste0(susie_output_dir, "/susie_", phenotype_name, "_", susie_prior, "_prior_sumstats.rds"))
  SusieResult_combined <- rbind(SusieResult_combined, SusieResult)
}

phenotype_name <- "WGBS_Flu"
susie_output_dir <- paste0(susie_dir, "/output/", phenotype_name, "/susie_", susie_prior, "_prior")
dir.create(susie_output_dir, showWarnings = F, recursive = T)
saveRDS(SusieResult_combined, paste0(susie_output_dir, "/susie_", phenotype_name, "_", susie_prior, "_prior_sumstats.rds"))


cat("Combine susie results for WGBS NI ... \n")
SusieResult_combined <- data.frame()
phenotype_list <-  paste0("WGBS_NI_full_E", 3:7)

for( phenotype_name in phenotype_list) {

  ## Load susie input data
  cat("Load", phenotype_name, "results ... \n")
  ## checking genes missing susie result
  susie_output_dir <- paste0(susie_dir, "/output/", phenotype_name, "/susie_", susie_prior, "_prior")
  SusieResult <- readRDS(paste0(susie_output_dir, "/susie_", phenotype_name, "_", susie_prior, "_prior_sumstats.rds"))
  SusieResult_combined <- rbind(SusieResult_combined, SusieResult)
}

phenotype_name <- "WGBS_NI"
susie_output_dir <- paste0(susie_dir, "/output/", phenotype_name, "/susie_", susie_prior, "_prior")
dir.create(susie_output_dir, showWarnings = F, recursive = T)
saveRDS(SusieResult_combined, paste0(susie_output_dir, "/susie_", phenotype_name, "_", susie_prior, "_prior_sumstats.rds"))
