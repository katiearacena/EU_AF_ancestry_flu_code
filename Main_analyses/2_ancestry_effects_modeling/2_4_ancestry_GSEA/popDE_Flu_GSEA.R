## Load libraries
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(fgsea)
library(msigdbr)
library(biomaRt)
library(org.Hs.eg.db)
library(dplyr)
library(ChIPseeker)
library(clusterProfiler)
library(stringi)
library(stringr)
library(reshape2)
library(data.table)
library(bsseq)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene


#############################################
## CREATE DIRECTORY STRUCTURE & LOAD DATA  ##
#############################################

## Set directory structure
folder = "2_ancestry_effects_modeling"
## Set directory 3 steps above script.
setwd('../../../')

system(paste0("mkdir -p Outputs/",folder,"/GSEA/"))
out_dir <- paste0("Outputs/",folder,"/GSEA/")


###############################
#### Get pathways for GSEA ####
###############################

## Get human pathways
m_df = msigdbr(species = "Homo sapiens")
## Check the available collections and sub-collections
m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
## Retrieve human genes for the hallmark collection gene sets
m_df = msigdbr(species = "Homo sapiens", category = "H")
human.path.list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

## Loop across all datasets
datatypes <- c("RNAseq", "ATACseq", "H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "WGBS")

for (i in 1:length(datatypes)){
    data_type_i <- datatypes[i]
    ## Print for reference
    print(data_type_i)

    system(paste0("mkdir -p ", out_dir, data_type_i, "/popDE_Flu/"))
    results_dir <- paste0(out_dir, data_type_i, "/popDE_Flu/")
    
    ## Set file and working directory (on local computer)
    file_name <- paste0("popDE_FLU.txt")
    working_dir <- paste0("Outputs/", folder,"/PopDE_analysis/", data_type_i, "/results/")

  
    ############################
    #### Load results table ####
    ############################
  
    ## for RNAseq data
    if(data_type_i =="RNAseq"){
        ## (because do not need to connect peaks to genes for RNAseq)
        ## Load DE results table.
        resultsALL <- read.csv(paste0(working_dir, file_name), sep=" ")
        ## get gene names from ensembe id (can also use biomaRt to do this but requires internet connection)
        gene_names <- read.table("Inputs/ref/ensemble_id_positions_with_gene_names.txt")
        gene_names <- select(gene_names, "gene_name", "Gene_ID")
        colnames(gene_names) <- c("gene_id", "symbol")
        resultsALL$symbol <-rownames(resultsALL)
        results_with_genes <- left_join(x= resultsALL, gene_names, by="symbol")

        ## Make gene rowname a column
        rownames(results_with_genes) <- results_with_genes$symbol
        results_with_genes$symbol <- results_with_genes$gene_id
        results_with_genes$gene_id <- NULL

        write.table(results_with_genes, file=paste0(results_dir, "resultsALL_with_genes.txt"))
  
    ## for WGBS data
    } else if (data_type_i == "WGBS"){

        ## Load file
        resultsALL <- read.csv(paste0(working_dir, file_name), sep=" ")
        ## Make peak id rowname a column
        resultsALL$peak_ids <- rownames(resultsALL)
    
        ## load granges object
        load("Inputs/counts_matrices/WGBS_filtered.counts.BSobj.RData")
        gr_full <- granges(BSobj.fit.nolow)
        ## filter to only include those  in final results file
        gr <- gr_full[gr_full@ranges@NAMES %in% rownames(resultsALL), ]
        ## add "chr" to gr file
        gr@seqnames <- sub("^", "chr", gr@seqnames)
  
    ##################################
    ### Connect CpG site with gene ###
    ##################################
    
        peakAnno <- annotatePeak(gr, tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")
    
        ## Save some cool plots
        pdf(paste0(results_dir,"/peakAnno_pichart_genomic_annotation.pdf"))
        plotAnnoPie(peakAnno)
        dev.off()
        
        pdf(paste0(results_dir,"peakAnno_barchart.pdf"))
        plotDistToTSS(peakAnno,
                    title="Distribution of transcription factor-binding loci\nrelative to TSS")
        dev.off()
        
        ## Make peakAnno file a dataframe
        peakAnno_df <- as.data.frame(peakAnno)

        ## save plots for only significant cpg sites
        resultsALL$peak_ids <- rownames(resultsALL)
        significant_hits <- resultsALL[resultsALL$qvals_Flu < .10, ]
        gr_sig <- gr_full[gr_full@ranges@NAMES %in% rownames(significant_hits), ]
        ## add "chr" to gr file
        gr_sig@seqnames <- sub("^", "chr", gr_sig@seqnames)
  
        significant_hits_peakAnno <- annotatePeak(gr_sig, tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")
        ## Save plots
        pdf(paste0(results_dir,"/peakAnno_pichart_genomic_annotation_SIG.pdf"))
        plotAnnoPie(significant_hits_peakAnno)
        dev.off()
        
        pdf(paste0(results_dir,"peakAnno_barchart_SIG.pdf"))
        plotDistToTSS(significant_hits_peakAnno,
                    title="Distribution of transcription factor-binding loci\nrelative to TSS (Significant PopDE Hits (q < .10)")
        dev.off()

        significant_hits_peakAnno_df <- as.data.frame(significant_hits_peakAnno)
        
        #############################################
        ## Merge gene symbol and id with cpg sites ##
        #############################################
        
        ## Extract relevant columns from peakAnno_df
        select_peakAnno_df <- data.frame(peakAnno_df$seqnames, peakAnno_df$start, peakAnno_df$end, peakAnno_df$SYMBOL, peakAnno_df$geneId)
        recomb_peaks <- data.frame(paste0(select_peakAnno_df$peakAnno_df.seqnames, "_", select_peakAnno_df$peakAnno_df.start), select_peakAnno_df$peakAnno_df.SYMBOL, select_peakAnno_df$peakAnno_df.geneId)
        colnames(recomb_peaks) <- c("peak_ids", "symbol", "gene_ids")
        recomb_peaks$peak_ids <- str_remove(recomb_peaks$peak_ids, "chr")
        
        ## Merge with results df
        results_with_genes <- full_join(x = resultsALL, y = recomb_peaks, by= "peak_ids")
        rownames(results_with_genes) <- results_with_genes$peak_ids
        write.table(results_with_genes, file=paste0(results_dir, "resultsALL_with_genes.txt"))

        ## write peak anno files
        write.table(significant_hits_peakAnno_df, file=paste0(results_dir, "significant_hits_peakAnno.txt"))
        write.table(peakAnno_df, file=paste0(results_dir, "peakAnno.txt"))
        
  ## for peak datatypes
    }else{
        ## Load file
        resultsALL <- read.csv(paste0(working_dir, file_name), sep=" ")
        ## Make peak id rowname a column
        resultsALL$peak_ids <- rownames(resultsALL)
        x <- as.data.frame(resultsALL$peak_ids)
        
        ###############################
        ### Connect peaks with gene ###
        ###############################
        
        str1 <- stri_replace_last_fixed(x$`resultsALL$peak_ids`, "_", "\t")
        str2 <- stri_replace_last_fixed(str1, "_", "\t")
        end <- colsplit(str2, '\t', names =  c('chr','start', 'end'))
        df <- as.data.frame(end)
        df_final <- print(df, quote=FALSE)
        
        ## Get GRanges object
        gr <- makeGRangesFromDataFrame(df_final)
        
        peakAnno <- annotatePeak(gr, tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")
        
        ## Save some cool plots
        pdf(paste0(results_dir, "peakAnno_pichart_genomic_annotation.pdf"))
        plotAnnoPie(peakAnno)
        dev.off()
        
        pdf(paste0(results_dir, "peakAnno_barchart.pdf"))
        plotDistToTSS(peakAnno,
                    title="Distribution of transcription factor-binding loci\nrelative to TSS")
        dev.off()
        
        ## Make peakAnno file a dataframe
        peakAnno_df <- as.data.frame(peakAnno)

        ## save plots for only significant cpg sites
        resultsALL$peak_ids <- rownames(resultsALL)
        significant_hits <- resultsALL[resultsALL$qvals_Flu < .10, ]
        significant_hits_x <- as.data.frame(significant_hits$peak_ids)
        
        ### Connect peaks with gene ###
        str1 <- stri_replace_last_fixed(significant_hits_x$`significant_hits$peak_ids`, "_", "\t")
        str2 <- stri_replace_last_fixed(str1, "_", "\t")
        end <- colsplit(str2, '\t', names =  c('chr','start', 'end'))
        significant_hits_df <- as.data.frame(end)
        significant_hits_df_final <- print(significant_hits_df, quote=FALSE)

        ## Get GRanges object
        significant_hits_gr <- makeGRangesFromDataFrame(significant_hits_df_final)
        significant_hits_peakAnno <- annotatePeak(significant_hits_gr, tssRegion=c(-3000, 3000),
                                TxDb=txdb, annoDb="org.Hs.eg.db")
        
        ## Save plots
        pdf(paste0(results_dir,"/peakAnno_pichart_genomic_annotation_SIG.pdf"))
        plotAnnoPie(significant_hits_peakAnno)
        dev.off()
        
        pdf(paste0(results_dir,"peakAnno_barchart_SIG.pdf"))
        plotDistToTSS(significant_hits_peakAnno,
                    title="Distribution of transcription factor-binding loci\nrelative to TSS (Significant PopDE Hits (q < .10)")
        dev.off()

        significant_hits_peakAnno_df <- as.data.frame(significant_hits_peakAnno)
        
        ############################################
        ## Merge gene symbol and id with peak_ids ##
        ############################################
        
        ## Extract relevant columns from peakAnno_df
        select_peakAnno_df <- data.frame(peakAnno_df$seqnames, peakAnno_df$start, peakAnno_df$end, peakAnno_df$SYMBOL, peakAnno_df$geneId)
        recomb_peaks <- data.frame(paste0(select_peakAnno_df$peakAnno_df.seqnames, "_", select_peakAnno_df$peakAnno_df.start, "_", select_peakAnno_df$peakAnno_df.end), select_peakAnno_df$peakAnno_df.SYMBOL, select_peakAnno_df$peakAnno_df.geneId)
        colnames(recomb_peaks) <- c("peak_ids", "symbol", "gene_ids")
        
        ## Merge with results df
        results_with_genes <- full_join(x = resultsALL, y = recomb_peaks, by="peak_ids")
        write.table(results_with_genes, file=paste0(results_dir, "resultsALL_with_genes.txt"))

        ## write peak anno files
        write.table(significant_hits_peakAnno_df, file=paste0(results_dir, "significant_hits_peakAnno.txt"))
        write.table(peakAnno_df, file=paste0(results_dir, "peakAnno.txt"))
    }
    
    ###############################
    ###### Format input file ######
    ###############################
    
    if(data_type_i =="WGBS"){
        ## Extract gene symbol and stat 
        ranks <- results_with_genes[, c("symbol", "stat")]
        ## Order t statistic in descending order (from + to - t statistic)
        ranks_t <- ranks[order(-(ranks$stat)), ]
        ranks_for_gsea <- setNames(ranks_t$stat, ranks_t$symbol)
    }else{

        ## Extract gene symbol and t statistic
        ranks <- results_with_genes[, c("symbol", "t")]
        ## Order t statistic in descending order (from + to - t statistic)
        ranks_t <- ranks[order(-(ranks$t)), ]
        ranks_for_gsea <- setNames(ranks_t$t, ranks_t$symbol)
    }
    
    ## File to save printed outputs to
    sink(paste0(results_dir, "filtering_stats.txt"))
    
    ## Total number of features
    print("Print total number of features before filtering")
    print(length(ranks_for_gsea))
    
    ## Number of unique genes
    print("Print number of unique genes")
    print(length(unique(names(ranks_for_gsea))))
    
    ## Number of duplicates before filtering
    print("Print the number of duplicates before filtering")
    print(anyDuplicated(names(ranks_for_gsea)))
    
    ## Save output
    sink()
    
    ## Remove duplicates by selecting the largest t-statistic value (the strongest association) per gene
    ranks_for_gsea_no_dups <- ranks_for_gsea[!duplicated(names(ranks_for_gsea))]
    ## Make sure this df is ordered
    ranks_for_gsea_no_dups <- ranks_for_gsea_no_dups[order(-(ranks_for_gsea_no_dups))]
    
    ## Number of duplicates after filtering
    anyDuplicated(names(ranks_for_gsea_no_dups))
    ## Number of features left after filtering
    length(ranks_for_gsea_no_dups)

    
    ###############################
    ######## Perform GSEA ########
    ###############################
    
    fgseaRes <- fgsea(pathways = human.path.list,
                        stats = ranks_for_gsea_no_dups, 
                        minSize=15, 
                        maxSize=500, 
                        nperm=100000)
    
    #Save results.
    fwrite(fgseaRes, file=paste0(results_dir, "fgsea_results.txt"), sep="\t", sep2=c("", " ", ""))
    
}
