
## pvalue to zscores
p2z <- function(pval, beta=NULL, two.sided=TRUE, limit=.Machine$double.xmin) {
  if( !is.null( limit ) ){
    pval[which(pval < limit )] <- limit ## set lower limit to avoid Inf/-Inf zscores
  }
  if(two.sided == FALSE){
    z <- qnorm(pval, lower.tail = FALSE)
  } else{
    z <- qnorm(pval/2, lower.tail = FALSE)
  }
  if ( !is.null( beta) ) {
    z <-  z * sign( beta )
  }
  return(z)
}

## Get genotype dosage from vcf files
get_genotype_dosage <- function(vcf_file, gene_bed, window_size = 1e5,
                                snpID_list, snpPOS_list,
                                path_vcftools = "vcftools", filterMAF = FALSE,
                                tmpdir = "/scratch/midway2/kaixuan/tmp/"){
  library(vcfR)
  library(data.table)

  ## Helper function to convert genotype into dosage.
  .genoDosage <- function(x){
    x <- as.character(x)
    x[grep("\\.", x)] <- NA
    return( as.numeric(stringr::str_count(x,"1")) )
  }

  tmpVcf <- tempfile(tmpdir = tmpdir, fileext = ".vcf")

  if(!missing(snpID_list)){
    # use SNP IDs
    cat("Load genotype using SNP IDs ...\n")
    tmp_snpID_list <- tempfile(tmpdir = tmpdir, fileext = ".snpIDs.txt")
    fwrite(as.data.frame(snpID_list), tmp_snpID_list, col.names = F)
    vcftools_cmd <- paste(path_vcftools, "--gzvcf", vcf_file, "--snps", tmp_snpID_list, "--recode --out", tmpVcf)
    cat(vcftools_cmd, '\n')
    system(vcftools_cmd)
    unlink(tmp_snpID_list) # remove the temp file to free space.

  }else if(!missing(snpPOS_list)){
    # use positions
    cat("Load genotype using SNP positions ...\n")
    tmp_snpPOS_list <- tempfile(tmpdir = tmpdir, fileext = ".snpPOS.txt")
    snpPOS_list <- as.data.frame(snpPOS_list)
    colnames(snpPOS_list) <- c("CHROM", "POS")
    fwrite(snpPOS_list, tmp_snpPOS_list, col.names = T, sep = "\t")
    vcftools_cmd <- paste(path_vcftools, "--gzvcf", vcf_file, "--positions", tmp_snpPOS_list, "--recode --out", tmpVcf)
    cat(vcftools_cmd, '\n')
    system(vcftools_cmd)
    unlink(tmp_snpPOS_list) # remove the temp file to free space.

  }else{
    # use positions
    cat("Load genotype using cis window ...\n")
    SNP_range <- data.frame(chr = gene_bed$chr,
                            start = max(0, as.integer(gene_bed$start) - window_size),
                            end = as.integer(gene_bed$end) + window_size)
    vcftools_cmd <- paste(path_vcftools, "--gzvcf", vcf_file, "--chr", SNP_range$chr, "--from-bp", SNP_range$start, "--to-bp", SNP_range$end, "--recode --out", tmpVcf)
    cat(vcftools_cmd, '\n')
    system(vcftools_cmd)
    # system(paste0("zcat ",vcf_file," | awk 'NR==1 {print $0} (/^#CHROM/){print $0} (!/^#/ && $1 == ", SNP_range$chr,
    #               " && $2 >= ", SNP_range$start," && $2 <= ",SNP_range$end," ) {print $0}' > ",tmpVcf, ".recode.vcf"))

  }

  # Extract genotype using vcfR
  cat("Extract genotype using vcfR ...\n")
  geno.vcf <- try(read.vcfR(file = paste0(tmpVcf, ".recode.vcf"), verbose = F), silent = T)
  unlink(paste0(tmpVcf, ".recode.vcf")) # remove the temp file to free space.

  ## check if genotype available for this peak
  if( nrow(geno.vcf@fix) == 0 ){
    cat(paste("Warning: no SNPs found! \n"))
    return(NULL)
  }

  ## Get genotype as dosage format and filter for MAF
  if( geno.vcf@gt[1,"FORMAT"] == "DS" ){
    ## Directly extract dosage
    geno <- if(nrow(geno.vcf@fix) == 1){
      t(apply(extract.gt(geno.vcf, element = 'DS' ),2,as.numeric ) )
    }else{
      apply(extract.gt(geno.vcf, element = 'DS' ),2,as.numeric )
    }
    rownames(geno) <- geno.vcf@fix[,"ID"]
    # rownames(geno) <- paste0(geno.vcf@fix[,"CHROM"], ":", geno.vcf@fix[,"POS"])

    ## filter out any genotype that has MAF<0.05
    if(filterMAF){
      MAF <- apply(geno,1,function(x) !any(table(round(x) ) > 0.95*ncol(geno)) )
      geno <- geno[MAF,]
      geno.vcf <- geno.vcf[MAF,]
    }
  }else if( geno.vcf@gt[1,"FORMAT"] == "GT" ){
    geno.vcf <- geno.vcf[is.biallelic(geno.vcf),]
    ## get genotype as Dosage
    tmp_geno <- extract.gt(geno.vcf, element = 'GT' )
    geno <- t( apply( tmp_geno, 1, .genoDosage ) )
    colnames(geno) <- colnames(tmp_geno)
    rownames(geno) <- geno.vcf@fix[,"ID"]
    # rownames(geno) <- paste0(geno.vcf@fix[,"CHROM"], ":", geno.vcf@fix[,"POS"])

    ## filter out any genotype that has MAF<0.05
    if(filterMAF){
      MAF <- apply(geno,1,function(x) !any(table(x) > 0.95*ncol(geno)) )
      geno <- geno[MAF,]
      geno.vcf <- geno.vcf[MAF,]
    }
  }

  if( nrow(geno) == 0 ){ return(NULL) } ## skip this iteration if no genotype left after filtering
  return(list(geno = geno, geno.vcf = geno.vcf))

}


## Extract genotype dosage from vcf file
get_genotype_dosage_simple <- function(vcf_file, gene_bed, window_size = 1e5, snpID_list, path_vcftools = "vcftools"){

  ## Helper function to convert genotype into dosage.
  .genoDosage <- function(x){
    x <- as.character(x)
    x[grep("\\.", x)] <- NA
    return( stringr::str_count(x,"1") )
  }

  tmpVcf <- tempfile(fileext = ".vcf")

  if(!missing(snpID_list)){
    # use SNP IDs
    cat("load genotype using SNP IDs ...\n")
    tmp_snpID_list <- tempfile(fileext = ".snpIDs.txt")
    fwrite(as.data.frame(snpID_list), tmp_snpID_list, col.names = F)
    system(paste(path_vcftools, "--gzvcf", vcf_file, "--snps", tmp_snpID_list, "--recode --out", tmpVcf))
    unlink(tmp_snpID_list) # remove the temp file to free space.
  }else{
    # use positions
    cat("load genotype using cis window ...\n")
    SNP_range <- data.frame(chr = gene_bed$chr,
                            start = max(0, as.integer(gene_bed$start) - window_size),
                            end = as.integer(gene_bed$end) + window_size)
    system(paste(path_vcftools, "--gzvcf", vcf_file, "--chr", SNP_range$chr, "--from-bp", SNP_range$start, "--to-bp", SNP_range$end, "--recode --out", tmpVcf))
  }

  gene_vcf <- as.data.frame(fread(paste0(tmpVcf, ".recode.vcf")))

  colnames(gene_vcf) <- gsub("#", "", colnames(gene_vcf))
  gene_vcf <- gene_vcf[gene_vcf$ID != ".", ]

  geno_info <- gene_vcf[,c(1:5)]
  tmp_geno <- gene_vcf[,-c(1:9)]
  geno <- t( apply( tmp_geno ,1, .genoDosage ) )
  colnames(geno) <- colnames(tmp_geno)
  # rownames(geno) <- paste0(gene_vcf[,"CHROM"], ":", gene_vcf[,"POS"])
  rownames(geno) <- gene_vcf$ID

  unlink(paste0(tmpVcf, ".recode.vcf")) # remove the temp file to free space.

  return(list(geno = geno, geno_info = geno_info))
}


prepare_susie_input_data <- function(vcf_file, phenotype_matrix, covariates_matrix, QTL_res, GENE, outdir) {

  # peak_bed <- phenotype[phenotype$gene == GENE, 1:4]

  gene_QTL_res <- dplyr::filter(QTL_res, gene == GENE)
  gene_QTL_res$zval <- p2z(gene_QTL_res$pvalue, gene_QTL_res$beta)

  # phenotype vector: n x 1
  y <- phenotype_matrix[, GENE]

  # covariate matrix: n x (k+1), including intercept
  Z <- cbind(intercept = 1, covariates_matrix)

  # genotype matrix: n x p, p: number of SNPs
  geno.l <- get_genotype_dosage(vcf_file, gene_QTL_res$SNP, GENE, filename_log)
  X <- t(as.matrix(geno.l$geno))

  ## select common SNPs and sort SNPs
  SNP_order <- intersect(colnames(X), gene_QTL_res$SNP)
  X <- X[, match(SNP_order, colnames(X))]
  gene_QTL_res <- gene_QTL_res[match(SNP_order, gene_QTL_res$SNP), ]

  est <- gene_QTL_res[,c("SNP","beta","se")]
  zval <- gene_QTL_res[,c("SNP","zval")]

  ## sort individuals
  X <- X[names(y), ]
  Z <- Z[names(y), ]

  ######## Removing covariates from X ########
  A <- crossprod(Z)

  # chol decomposition for (Z'Z)^(-1)
  R <- chol(solve(A))
  W <- R %*% t(Z) %*% X

  # Remove Covariates from X
  Xhat <- as.matrix(X - Z %*% crossprod(R,W))

  # LD (R)
  R_LD <- cor(Xhat)

  ######## Removing covariates from Y ########
  A <- crossprod(Z)
  SZy <- solve(A, drop(y %*% Z))
  yhat <- y - drop(Z %*% SZy)

  ###############################################################
  ## Codes above are the implementation of the following formula:
  # yhat = y - Z*((Z'*Z)\(Z'*y))
  # Xhat = X - Z*((Z'*Z)\(Z'*X))
  ###############################################################

  dir.create(paste0(outdir, "/INPUT/", GENE), showWarnings = F, recursive = T)
  write.table(est, paste0(outdir, "/INPUT/", GENE, "/est.dat"), sep = "\t", col.names = F,row.names = F,quote = F)
  write.table(zval, paste0(outdir, "/INPUT/", GENE, "/zval.dat"), sep = "\t", col.names = F,row.names = F,quote = F)
  write.table(R_LD, paste0(outdir, "/INPUT/", GENE, "/R_LD.dat"), sep = "\t", col.names = F,row.names = F,quote = F)
  save(Xhat, yhat, file = paste0(outdir, "/INPUT/", GENE, "/individual_data.RData"))

}


## Combine Susie result of each gene into a full summary stats for all gene-SNP pairs
combine_susie_results <- function(fitted_susie_allGenes, geneIDs, num_cores = 5){
  registerDoParallel(num_cores)

  fitted_susie_allGenes[sapply(fitted_susie_allGenes, is.null)] <- NULL

  if(missing(geneIDs)){
    geneIDs <- names(fitted_susie_allGenes)
  }

  SusieResult <- foreach( iGene = geneIDs, .combine = rbind) %dopar% {

    SusieResult_gene <- fitted_susie_allGenes[[iGene]]

    Zscore <- SusieResult_gene$z_score
    SNP <- names(Zscore)

    fitted_susie <- SusieResult_gene$fitted

    tmpOut <- data.frame(SNP = SNP, GENE = iGene, PIP = fitted_susie$pip, CS = NA, Zscore = Zscore)
    idx_inCS <- unlist(fitted_susie$sets$cs)
    if( !is.null(idx_inCS) ){
      tmpOut$CS[idx_inCS] <- rep(names(fitted_susie$sets$cs), lapply(fitted_susie$sets$cs, length))
    }

    tmpOut[order(tmpOut$PIP,decreasing = T),]

  }

  SusieResult$CS[is.na(SusieResult$CS)] <- 0
  SusieResult$CS <- as.integer(gsub("L", "", SusieResult$CS))

  stopImplicitCluster()

  return(SusieResult)
}


fit_susie_uniform_prior <- function(X, y, covariates, L = 10, plot = FALSE) {

  ## helper function to remove covariate effects
  .remove.covariate.effects <- function (X, Z, y) {
    # include the intercept term
    if (any(Z[,1]!=1)) Z = cbind(1, Z)
    A   <- forceSymmetric(crossprod(Z))
    SZy <- as.vector(solve(A,c(y %*% Z)))
    SZX <- as.matrix(solve(A,t(Z) %*% X))
    y <- y - c(Z %*% SZy)
    X <- X - Z %*% SZX
    return(list(X = X,y = y,SZy = SZy,SZX = SZX))
  }

  X <- as.matrix(X)

  ## remove the individuals with NA, and keep the rest of samples for this gene
  if(anyNA(y)){
    cat('Remove', length(which(is.na(y))), 'samples with missing values. \n')
    idv_nonNA <- names(y)[which(!is.na(y))]
    y <- y[idv_nonNA]
    X <- X[idv_nonNA, ]
    covariates <- covariates[idv_nonNA,]
  }

  snps_withNAs <- apply(X, 2, function(x) anyNA(x))
  if(any(snps_withNAs)){
    cat('Remove', length(snps_withNAs), 'variants with missing values. \n')
    X <- X[, !snps_withNAs]
  }

  if(!is.null(covariates)){
    # covariate matrix: n x (k+1), including intercept
    Z <- cbind(intercept = 1, covariates)
    out <- .remove.covariate.effects(X, Z, y)
    Xhat <- out$X
    yhat <- out$y
  }else{
    Xhat <- as.matrix(X)
    yhat <- y
  }

  ## compute PVE
  reg <- susieR::univariate_regression(Xhat, yhat)
  z_score <- reg$betahat/reg$sebetahat
  beta <- reg$betahat
  names(beta) <- colnames(Xhat)
  names(z_score) <- colnames(Xhat)
  top_idx <- which.max(abs(z_score))
  PVE <- var(Xhat[,top_idx] * reg$betahat[top_idx]) / var(yhat)
  cat("PVE =",PVE, "\n")


  ## fit susie with uniform prior
  # cat("Fit susie with uniform prior...\n")
  fitted <- susieR::susie(Xhat,
                          yhat,
                          L = L,
                          prior_weights = NULL,
                          estimate_residual_variance = TRUE,
                          estimate_prior_variance = TRUE)

  if(plot) {
    susie_plot(fitted, y="PIP")
  }

  return(list(fitted = fitted,
              z_score = z_score,
              PVE = PVE))

}

## Fit susie with individual level data
fit_susie_idv <- function(X, y, covariates=NULL, prior=NULL, L = 10, plot = FALSE) {

  ## helper function to remove covariate effects
  .remove.covariate.effects <- function (X, Z, y) {
    # include the intercept term
    if (any(Z[,1]!=1)) Z = cbind(1, Z)
    A   <- forceSymmetric(crossprod(Z))
    SZy <- as.vector(solve(A,c(y %*% Z)))
    SZX <- as.matrix(solve(A,t(Z) %*% X))
    y <- y - c(Z %*% SZy)
    X <- X - Z %*% SZX
    return(list(X = X,y = y,SZy = SZy,SZX = SZX))
  }

  X <- as.matrix(X)

  ## remove the individuals with NA, and keep the rest of samples for this gene
  if(anyNA(y)){
    # cat('Remove', length(which(is.na(y))), 'samples with missing values. \n')
    idv_nonNA <- names(y)[which(!is.na(y))]
    y <- y[idv_nonNA]
    X <- X[idv_nonNA,,drop=FALSE]
    covariates <- covariates[idv_nonNA,,drop=FALSE]
  }

  snps_withNAs <- apply(X, 2, function(x) anyNA(x))
  if(any(snps_withNAs)){
    # cat('Remove', length(snps_withNAs), 'variants with missing values. \n')
    X <- X[, !snps_withNAs, drop=FALSE]
  }

  if(!is.null(covariates)){
    # covariate matrix: n x (k+1), including intercept
    Z <- cbind(intercept = 1, covariates)
    out <- .remove.covariate.effects(X, Z, y)
    Xhat <- out$X
    yhat <- out$y
  }else{
    Xhat <- as.matrix(X)
    yhat <- y
  }

  ## compute PVE
  reg <- susieR::univariate_regression(Xhat, yhat)
  z_score <- reg$betahat/reg$sebetahat
  beta <- reg$betahat
  names(beta) <- colnames(Xhat)
  names(z_score) <- colnames(Xhat)
  top_idx <- which.max(abs(z_score))
  PVE <- var(Xhat[,top_idx] * reg$betahat[top_idx]) / var(yhat)
  # cat("PVE =",PVE, "\n")

  if(is.null(prior)){
    ## fit susie with uniform prior
    # cat("Fit susie with uniform prior...\n")
    fitted <- susieR::susie(Xhat,
                            yhat,
                            L = L,
                            prior_weights = NULL,
                            estimate_residual_variance = TRUE,
                            estimate_prior_variance = TRUE)
  }else{
    ## fit susie with prior
    # cat("Fit susie with prior...\n")
    fitted <- susieR::susie(Xhat,
                            yhat,
                            L = L,
                            prior_weights = prior,
                            estimate_residual_variance = TRUE,
                            estimate_prior_variance = TRUE)
  }

  if(plot) {
    susie_plot(fitted, y="PIP")
  }

  return(list(fitted = fitted,
              z_score = z_score,
              PVE = PVE))

}
