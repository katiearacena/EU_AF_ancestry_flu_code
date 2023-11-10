## Written on 1/12/23  by KA 🐣

## on cluster is named pt3_calc_overlap.R
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

##set working directory.
inDir="/project/lbarreiro/USERS/katie/EU_AF_ancestry_flu/local_ancestry/popDE_inputs/"
LA_genos <- data.frame(fread(paste0(inDir, "LA_matrix_for_popDE.txt")))
LA_genos$V1 <- NULL
LA_genos$LA_region <- paste0("chr", LA_genos$chr, "_",LA_genos$spos, "_", LA_genos$epos)
## remove old cols
LA_genos <- select(LA_genos, -c("chr", "spos", "epos"))

## Load results
datatypes <- c("RNAseq","ATACseq","H3K27ac","H3K27me3","H3K4me1", "H3K4me3", "WGBS")

for (i in 1:length(datatypes)){
  
  data_type_i <- datatypes[i]
  ## Print for reference
  print(data_type_i)
  
  ## load file
  regions <- data.frame(fread(paste0(inDir, data_type_i, "_intersected_LA_coordinates.bed")))
  ## re-format
  regions$extended_region <- paste0("chr", regions$V1, "_",regions$V2, "_", regions$V3)
  regions$feature <- paste0("chr", regions$V4, "_",regions$V5, "_", regions$V6)
  regions$LA_region <- paste0("chr", regions$V7, "_",regions$V8, "_", regions$V9)
  regions$bp_overlap <- paste0(regions$V10)
  
  coords <- select(regions, extended_region, feature, LA_region, bp_overlap)
  coords_g <- coords %>% group_by(extended_region, feature) %>% summarise_all(toString)
  
  ##make matrix built row by row
  ## set colnames
  new_matrix <- matrix(, nrow = nrow(coords_g), ncol = 37)
  colnames(new_matrix) <- c("feature", colnames(LA_genos))
  
  for (j in 1:nrow(coords_g)){
    

    if (j%%1000 == 0) {print(j)}
    
    if(max(str_count(coords_g[j, ]$LA_region, ','))==0){
      
      ## match LA region with genos.
      ## it will be 100% the genotype of this genotype since this is the only overlap
      new_matrix[j,2:37 ]  <- as.vector(as.matrix(LA_genos[LA_genos$LA_region == coords_g[j, ]$LA_region, ]))
      new_matrix[j,2:36 ] <- recode(new_matrix[j,2:36 ],`0` =  0,  `1` = .5, `2` = 1)
      ## write the original gene/peak  as well
      new_matrix[j, 1] <-coords_g[j, ]$feature 
      
    }else if(max(str_count(coords_g[j, ]$LA_region, ','))>0){
      
      #print("multiple overlaps")
      
      ## for each feature, split by comma
      overlapping_regions <- unlist(strsplit(coords_g[j, ]$LA_region, ","))
      ## remove any extra spaces
      overlapping_regions <- str_replace(overlapping_regions, " ", "")
      
      ## also get distance for each feature
      overlapping_distances_for_region <- unlist(strsplit(coords_g[j, ]$bp_overlap, ","))
      
      for(f in 1:length(overlapping_regions)){
        overlapping_regions_f <- overlapping_regions[f]
        bp_overlap_f <- overlapping_distances_for_region[f]
        
        ## get the genos from matrix for each feature
        if(f==1){
          tmp_matrix <-  LA_genos[LA_genos$LA_region == overlapping_regions_f, ]
        } else {
          tmp_matrix <- rbind(tmp_matrix, LA_genos[LA_genos$LA_region == overlapping_regions_f, ])
        }
      }

      tmp_matrix[ , 37] <- as.numeric(overlapping_distances_for_region)
      colnames(tmp_matrix)[37] <- "num_bp_overlap"
      ## recoded alleles
      tmp_matrix[,1:35 ] <- recode(as.matrix(tmp_matrix[,1:35 ]),`0` =  0,  `1` = .5, `2` = 1)
      
      ## now need to calculate overall AF and EU ancestry for this region
      tmp_matrix$weight = tmp_matrix$num_bp_overlap/sum(tmp_matrix$num_bp_overlap)
      ## multiple genos by weights and take sum of rows to get 1 meta value
      weighted_matrix <- colSums(tmp_matrix[, 1:35] * tmp_matrix$weight)
      ## now write as full line to add to full matrix
      weighted_matrix <- t(as.data.frame(weighted_matrix))
      weighted_matrix <- as.data.frame(cbind(weighted_matrix, coords_g[j, ]$feature, coords_g[j, ]$LA_region))
      colnames(weighted_matrix)[36] <- "feature"
      colnames(weighted_matrix)[37] <- "LA_region"
      ## reorder
      weighted_matrix <- weighted_matrix %>% select(colnames(new_matrix))
      new_matrix[j, 1:37]  <- as.matrix(weighted_matrix)
      
    }
  }
  write.table(new_matrix, file=paste0(data_type_i, "_matrix_for_LA_v2.txt"))
  print(paste0(data_type_i, " completed"))
}
