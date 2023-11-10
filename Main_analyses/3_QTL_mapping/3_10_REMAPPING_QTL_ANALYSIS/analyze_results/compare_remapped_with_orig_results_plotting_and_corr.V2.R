##############################

## correlate BEST snps from remapped and original data
library(dplyr)
library(ggpubr)
library(data.table)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

args<-commandArgs(TRUE)
## Declare condition equal to "NI or "Flu" 
condition <- args[1]
## Declare data type.
datatype <- args[2]
## Declare number of expression PCs to regress out.
expPCs <- args[3]

if(length(args)==0)
{
  print("WARNING: No arguments supplied.")
}
print(args)


## can make these args
#condition <- "NI"
#datatype <- "RNAseq"
#expPCs <- "4"

dir <- "/project2/lbarreiro/users/katie/EU_AF_ancestry_flu/"
out_dir <- "/project2/lbarreiro/users/katie/EU_AF_ancestry_flu/paper_pipeline/Outputs/3_QTL_mapping/SNP-QTL_mapping-REMAPPED_cts/"
  
## load data
if (datatype=="ATACseq"){
  remapped_best <-  read.table(paste0(dir, "/paper_pipeline/Outputs/3_QTL_mapping/SNP-QTL_mapping-REMAPPED_cts/", datatype, "/", condition, "/", datatype,"_", condition, "_", expPCs,"/results_best_SNPs_with_qval.txt"))
} else {
  remapped_best <-  read.table(paste0(dir, "/paper_pipeline/Outputs/3_QTL_mapping/SNP-QTL_mapping-REMAPPED_cts/", datatype, "/", condition, "/", datatype,"_", condition, "/results_best_SNPs_with_qval.txt"))
}
orig_best <- read.table(paste0(dir, "paper_pipeline/Outputs/3_QTL_mapping/SNP-QTL_mapping/", datatype, "/", condition, "/", datatype, "_", condition, "/results_best_SNPs_with_qval.txt"))
 
## subset
remapped_best_sig <- remapped_best[remapped_best$qvalues <.10, ]$gene
remapped_best_betas <- select(remapped_best, gene, beta)

orig_best_sig <- orig_best[orig_best$qvalues <.10, ]$gene
orig_best_betas <- select(orig_best, gene, beta)

t1 <- print(paste0("#sig hits in remapped (q<.10)=", length(remapped_best_sig)))
t2 <- print(paste0("#sig hits in original (q<.10)=", length(orig_best_sig)))
t3 <- print(paste0("#sig hits shared=", length(intersect(remapped_best_sig, orig_best_sig))))

## compare effect sizes of best snps
df <- merge(remapped_best_betas, orig_best_betas, by.x="gene", by.y="gene")

p1<- ggscatter(df, x="beta.x", y="beta.y",add = "reg.line",  add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) + 
  stat_cor(method = "pearson") + ylab("original beta") +xlab("remapped beta") +ggtitle("BEST peaks, all peaks")

## subsetting on those that are QTL in EITHER
sig_features_in_either <- union(remapped_best_sig, orig_best_sig)

df_sig <- df[df$gene %in% sig_features_in_either, ]
p3 <- ggscatter(df_sig, x="beta.x", y="beta.y",add = "reg.line",  add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) + 
  stat_cor(method = "pearson") + ylab("original beta") +xlab("remapped beta")+ ggtitle("BEST peaks, q<.10 filter in EITHER")

best_text <- textGrob(paste0(datatype,", ",condition,", expression PCs used in both=",expPCs, "\n", t1, "\n", t2, "\n", t3, "\n"))

#png(filename=paste0("comparison_of_original_v_remapping_", datatype, condition, expPCs, ".png"), width=8, height=10, units="in", res=800)
#grid.arrange(best_text, all_text, p1, p1_all, p3, p_all, nrow=3, ncol=2)
#grid.arrange(best_text, best_text, p1, p1, p3, p3, nrow=3, ncol=2)
#dev.off()

###############################
### same but with all snps  ###
###############################

## need to do this part on cluster due to memory limit ****
## load data
if (datatype=="ATACseq"){
  remapped_all <-  as.data.frame(fread(paste0(dir, "paper_pipeline/Outputs/3_QTL_mapping/SNP-QTL_mapping-REMAPPED_cts/", datatype, "/", condition, "/", datatype,"_", condition, "_", expPCs,"/raw_results/result_original.txt")))
} else {
  remapped_all <-  as.data.frame(fread(paste0(dir, "paper_pipeline/Outputs/3_QTL_mapping/SNP-QTL_mapping-REMAPPED_cts/", datatype, "/", condition, "/", datatype,"_", condition, "/raw_results/result_original.txt")))
}

orig_all <- as.data.frame(fread(paste0(dir, "paper_pipeline/Outputs/3_QTL_mapping/SNP-QTL_mapping/", datatype, "/", condition, "/", datatype,"_", condition, "/raw_results/result_original.txt")))

## subset
## load appropriate q value equivalent pvalue threshold calculated from original mapping results
thresholds <-readRDS(paste0(dir, "Make_figures/QTL_integration_thresholds.q.10.RData"))
threshold <- thresholds[paste0(condition, "_", datatype, "_q.10")]

remapped_all_sig <- remapped_all[remapped_all$pvalue <threshold, ]$gene
remapped_all$snp_gene_pair <- paste0(remapped_all$snp, "_", remapped_all$gene)
remapped_all_betas <- select(remapped_all, snp_gene_pair, beta)
remapped_all_sig_to_count <- remapped_all[remapped_all$pvalue <threshold, ]$snp_gene_pair

orig_all_sig <- orig_all[orig_all$pvalue <threshold, ]$gene
orig_all$snp_gene_pair <- paste0(orig_all$snp, "_", orig_all$gene)
orig_all_betas <- select(orig_all, snp_gene_pair, beta)
orig_all_sig_to_count <- orig_all[orig_all$pvalue <threshold, ]$snp_gene_pair

t1.a <- print(paste0("#sig hits in remapped (q<.10 equiv. P)=", length(remapped_all_sig_to_count)))
t2.a <- print(paste0("#sig hits in original (q<.10 equiv. P)=", length(orig_all_sig_to_count)))
t3.a <- print(paste0("#sig hits shared=", length(intersect(remapped_all_sig_to_count, orig_all_sig_to_count))))

## compare effect sizes of best snps
df2 <- full_join(remapped_all_betas, orig_all_betas, by="snp_gene_pair")

#p1_all<- ggscatter(df2, x="beta.x", y="beta.y",add = "reg.line",  add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) + 
#  stat_cor(method = "pearson") + ylab("original beta") +xlab("remapped beta") +ggtitle("ALL peaks, all peaks")

## subsetting on those that are QTL in EITHER
sig_features_in_either2 <- union(remapped_all_sig_to_count, orig_all_sig_to_count)

df2_sig <- df2[df2$snp_gene_pair %in% sig_features_in_either2, ]
p3_all <- ggscatter(df2_sig, x="beta.x", y="beta.y",add = "reg.line",  add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) + 
  stat_cor(method = "pearson") + ylab("original beta") +xlab("remapped beta")+ ggtitle("ALL peaks, q<.10 equiv in EITHER")

all_text <- textGrob(paste0(datatype,", ",condition,", expression PCs used in both=",expPCs, "\n", t1.a, "\n", t2.a, "\n", t3.a, "\n"))

p1_vec <- cor.test(df2$beta.x, df2$beta.y,method = "pearson")
p1_all_corr_text <- textGrob(paste0("corr. coef = ", p1_vec$estimate, "\n p-value=", p1_vec$p.value))

p3_vec <- cor.test(df2_sig$beta.x, df2_sig$beta.y,method = "pearson")
p3_all_corr_text <- textGrob(paste0("corr. coef = ", p3_vec$estimate, "\n p-value=", p3_vec$p.value))

#png(filename=paste0(out_dir, "comparison_of_original_v_remapping_", datatype, condition, expPCs, ".png"), width=10, height=10, units="in", res=800)
#grid.arrange(best_text, all_text, p1, p1_all_corr_text, p3, p3_all, nrow=3, ncol=2)
#dev.off()

#######################################################################################
### plotting the effect size of the best snp AMONG 2 comparisons and plot that snp  ###
#######################################################################################

## select best snp for each feature

remapped_all_pt3 <- select(remapped_all, -c(V1, statistic, FDR))
orig_all_pt3 <- select(orig_all, -c(V1, statistic, FDR))

## join
df3 <- full_join(remapped_all_pt3, orig_all_pt3, by="snp_gene_pair")

event.sort_remapped_first<-df3[order(df3[,"gene.x"],df3[,"pvalue.x"]),]
event.sort_remapped_first_ordered <- event.sort_remapped_first[!duplicated(event.sort_remapped_first$gene.x),]

event.sort_orig_first<-df3[order(df3[,"gene.x"],df3[,"pvalue.y"]),]
event.sort_orig_first_ordered <- event.sort_remapped_first[!duplicated(event.sort_orig_first$gene.x),]

event.sort_remapped_first_ordered_PLOT<- ggscatter(event.sort_remapped_first_ordered, x="beta.x", y="beta.y",add = "reg.line",  add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) + 
  stat_cor(method = "pearson") + ylab("original beta") +xlab("remapped beta") +ggtitle("BEST snp among, sorted remapped pvals")

event.sort_orig_first_ordered_PLOT<- ggscatter(event.sort_orig_first_ordered, x="beta.x", y="beta.y",add = "reg.line",  add.params = list(color = "blue", fill = "lightgray"), conf.int = TRUE) + 
  stat_cor(method = "pearson") + ylab("original beta") +xlab("remapped beta") +ggtitle("BEST snp among, sorted orig pvals")

text_v2 <- textGrob(paste0(datatype,", ",condition,", expression PCs used in both=",expPCs, "\n plotting the best snp for each feature\n top: sorted remapped pvalues \n bottom: sorted original pvalues"))


png(filename=paste0(out_dir, "comparison_of_original_v_remapping_", datatype, condition, expPCs, "v2.png"), width=15, height=10, units="in", res=800)
grid.arrange(best_text, all_text, text_v2, p1, p1_all_corr_text, event.sort_remapped_first_ordered_PLOT, p3, p3_all, event.sort_orig_first_ordered_PLOT, nrow=3, ncol=3)
dev.off()
