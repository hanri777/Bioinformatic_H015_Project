# Adding Gene name columns to count_mat_df_orig to DESeqRes_df and best genes
# resource_df == DESeqRes_df or best_genes
# origin_count_matrix_df == count_mat_df_orig

func_create_gene_names_col <- function(resource_df, origin_count_matrix_df){
  

# Creating copy of resource_df (DESeqRes_df)
  #copy_best_genes <- best_genes
  copy_resource_df <- resource_df
  #copy_best_genes$Gene <- NA
  copy_resource_df$Gene <-NA
 #copy_count_mat_df_orig <- count_mat_df_orig
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

best_genes_names <- func_create_gene_names_col(best_genes, count_mat_df_orig)
best_genes_names

DESeqRes_df_Genenames <- func_create_gene_names_col(DESeqRes_df, count_mat_df_orig)
head(DESeqRes_df_Genenames)

#-------------

# Adding Gene name columns to count_mat_df_orig to DESeqRes_df and best genes
# resource_df == DESeqRes_df or best_genes
# origin_count_matrix_df == count_mat_df_orig

#func_create_gene_names_col <- function(resource_df, origin_count_matrix_df){

# Creating copy of resource_df (DESeqRes_df)
#  copy_resource_df <- resource_df
copy_best_genes <- best_genes
copy_best_genes$Gene <- NA
copy_count_mat_df_orig <- count_mat_df_orig

for(rows_best_genes in rownames(copy_best_genes)){
  for(rows_count_mat_df_orig in rownames(copy_count_mat_df_orig)){
    if(rows_best_genes == rows_count_mat_df_orig){
      copy_best_genes[rows_best_genes, "Gene"] = 
        copy_count_mat_df_orig[rows_count_mat_df_orig, "Gene"]
    }
  }
}

#  return(copy_resource_df)
#}

#best_genes_genes_names <- func_create_gene_names_col(best_genes, count_mat_df_orig)
#best_genes_genes_names

#DESeqRes_df_Genenames <- func_create_gene_names_col(DESeqRes_df, count_mat_df_orig)
#head(DESeqRes_df_Genenames)