DATA_TYPES=(RNAseq ATACseq H3K27ac H3K27me3 H3K4me1 H3K4me3 WGBS)
LEN=${#DATA_TYPES[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
       DAT=${DATA_TYPES[$NUM]}

        sed -i 's/_NI//g' ${DAT}_matrix_for_LA_v2.txt
        sed -i 's/_Flu//g' ${DAT}_matrix_for_LA_v2.txt
        sed -i 's/_Mock//g' ${DAT}_matrix_for_LA_v2.txt
        
    done

    
## add ensemble id to RNAseq

module load R
R

library(dplyr)

gene_pos <- read.table("/project/lbarreiro/USERS/katie/EU_AF_ancestry_flu/paper_pipeline/Inputs/QTL_mapping/RNAseq_positions.txt")
gene_pos$feature <- paste0(gene_pos$chromosome, "_", gene_pos$S1, "_", gene_pos$S2)

RNAseq_LA_matrix <- read.table("/project/lbarreiro/USERS/katie/EU_AF_ancestry_flu/local_ancestry/popDE_inputs/outputs/RNAseq_matrix_for_LA_v2.txt")

RNAseq_pos <- full_join(RNAseq_LA_matrix, gene_pos, by=c("feature" = "feature"))
RNAseq_pos <- select(RNAseq_pos, -c(chromosome, S1, S2))


write.table(RNAseq_pos, "RNAseq_matrix_for_LA_v2.txt")
