# Installing packages
install.packages("Rtsne")

# Loading packages
library(DESeq2)
library(vsn)
library(readxl)
library(ggplot2)
library(tidyverse)
library(dplyr) 
library(grid)
library(ComplexHeatmap)
library(ggrepel)
library(EnhancedVolcano)
library(pheatmap)
library(Rtsne)

# Loading and converting data to data frame
# Count Matrix
count_matrix <- read_xlsx('count_meta_origin_files/H015_Reference_Genome_Gu_counts.xlsx')
count_mat_df_orig <- as.data.frame(count_matrix)
class(count_mat_df_orig)
str(count_mat_df_orig)
dim(count_mat_df_orig)
head(count_mat_df_orig)

# Set Ens_ID column as row names
rownames(count_mat_df_orig) <- count_mat_df_orig$Ens_ID

# Deleting columns which are not samples and also not numeric
count_mat_df <- count_mat_df_orig[, 5 : length(count_mat_df_orig)]
dim(count_mat_df)
head(count_mat_df)
class(count_mat_df)

# Sorting columns
sorted_colnames <- sort(colnames(count_mat_df))
count_mat_df <- count_mat_df[, sorted_colnames, drop = FALSE]

# Meta Data
metadata <- read_xlsx('count_meta_origin_files/MetaShort_HIPO_15_negative_and_positive_samples.xlsx')
metadata_df <- as.data.frame(metadata)
class(metadata_df)

# Set Samples column as row names
rownames(metadata_df) <- metadata_df$Samples
class(metadata_df)
str(metadata_df)
dim(metadata_df)
head(metadata_df)

# Basic Quality Check of Feature Count Matrix and Meta Data
# Check if row names in metadata_df matches col names in count data
all(colnames(count_mat_df) %in% rownames(metadata_df))
# Checking order of row names and column names
all(rownames(metadata_df) == colnames(count_mat_df))

# Building DESeq Data Set
dds <- DESeqDataSetFromMatrix(countData = count_mat_df, 
                              colData = metadata_df,
                              design = ~ status)
dds

# Removing Low Counts Reads Genes 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

#Differential Analysis of Gene (compare treated with untreated)

# Set reference value for DEG analisys 
dds$status <- relevel(dds$status, ref = "negative")

# Perform Differential Analysis of Gene
deg <- DESeq(dds)

# Getting results and sorting as data frame
DESeqRes <- results(deg, alpha = 0.05)
DESeqRes <- DESeqRes[!is.na(DESeqRes$padj),]
DESeqRes <- DESeqRes[order(DESeqRes$padj),]
summary(DESeqRes)
DESeqRes_df <- as.data.frame(DESeqRes)
class(DESeqRes_df)

#Getting idea about best Genes
best_genes <- DESeqRes_df %>%
  arrange(padj) %>%
  head(20)
best_genes
class(best_genes)

#-------------------
# Creating function  
# Adding Gene name columns to count_mat_df_orig to DESeqRes_df and best genes
# resource_df == DESeqRes_df or best_genes
# origin_count_matrix_df == count_mat_df_orig
  
func_create_gene_names_col <- function(resource_df, origin_count_matrix_df){
    
  # Creating copy of resource_df (DESeqRes_df)
  copy_resource_df <- resource_df
  # Creating new empty column
  copy_resource_df$Gene <-NA
  # Creatubg copy of origin_count_matrix_df (count_mat_df_orig)
  copy_origin_count_matrix_df <- origin_count_matrix_df
    
  for(rows_resource_df in rownames(copy_resource_df)){
    for(rows_origin_count_matrix_df in rownames(copy_origin_count_matrix_df)){
      if(rows_resource_df == rows_origin_count_matrix_df){
        copy_resource_df[rows_resource_df, "Gene"] = 
          copy_origin_count_matrix_df[rows_origin_count_matrix_df, "Gene"]
      }
    }
  }
    
  return(copy_resource_df)
}

#------------
# Calling function  
best_genes_names <- func_create_gene_names_col(best_genes, count_mat_df_orig)
best_genes_names
  
DESeqRes_df_Genenames <- func_create_gene_names_col(DESeqRes_df, count_mat_df_orig)
head(DESeqRes_df_Genenames)  
  
#------------
# Writing to .csv file
write.csv(DESeqRes_df_Genenames, "output_files/DESeqResult_df_Genenames.csv")
write.csv(best_genes_names, "output_files/best_genes_names.csv")


#-------------------
### Plots

##Quality Check For RNA-Sec Data

# 1 PCA plot
vsd <- vst(deg, blind = FALSE)
plotPCA(vsd, intgroup = "status")

# 2 Size factor estimation
sizeFactors(deg)

# 3 Dispersion Plot
plotDispEsts(deg)

#MA plot
plotMA(DESeqRes) 

#------------------
# Volcono Plot
EnhancedVolcano(DESeqRes_df, 
                lab = rownames(DESeqRes_df), 
                x = 'log2FoldChange', 
                y = 'padj', 
                subtitle = 'Positive vs Negative', 
                labSize = 3, 
                pCutoff = 0.05,
                FCcutoff = 1,
                drawConnectors = TRUE,
                max.overlaps = 300)

#-----------------------
# Same Volcono Plot with gene names
EnhancedVolcano(DESeqRes_df, 
                lab = DESeqRes_df_Genenames$Gene, 
                x = 'log2FoldChange', 
                y = 'padj', 
                subtitle = 'Positive vs Negative', 
                labSize = 3, 
                pCutoff = 0.05,
                FCcutoff = 1,
                drawConnectors = TRUE,
                max.overlaps = 300)

# ----------------------
# Heat Map
top_genes <- DESeqRes_df %>% 
  arrange(padj) %>%
  head(30)

# Normalize deg data (DESec(dds))
mat <- counts(deg, normalized = T)[rownames(top_genes),]
head(mat, 5)

# Get z values
mat.z <- t(apply(mat, 1, scale))
head(mat.z, 5)

#Fix names of columns
colnames(mat.z) <- rownames(metadata_df)
head(mat.z, 5)

#Build Heatmap
Heatmap(mat.z, cluster_rows = T, cluster_columns = T, 
        column_labels = colnames(mat.z), 
        row_labels = top_genes$baseMean)

#-----------------------------
# Alternative better HeatMap
rltfd <- rlog(deg, blind=FALSE)
DEG.idx <- which(DESeqRes$padj <= 0.05 & 
                   abs(DESeqRes$log2FoldChange) > 1)
DESeqRes[DEG.idx,]

pheatmap(assay(rltfd)[DEG.idx,], 
         treeheight_row = 0, 
         treeheight_col = 0, 
         scale = "row")

# Alternative better HeatMap with gene names
DEG.idx <- which(DESeqRes_df_Genenames$padj <= 0.05 & 
                   abs(DESeqRes_df_Genenames$log2FoldChange) > 1)
DESeqRes_df_Genenames[DEG.idx,]

pheatmap(assay(rltfd)[DEG.idx,], 
         treeheight_row = 0, 
         treeheight_col = 0, 
         scale = "row",
         labels_row = DESeqRes_df_Genenames$Gene
         )


