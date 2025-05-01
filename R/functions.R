# remove genes with 80% zeroes and na rows
remove_sparse_rows <- function(df, threshold = 0.8) {
  df[rowMeans(df == 0, na.rm = TRUE) < threshold & !apply(df, 1, function(x) all(is.na(x))), , drop = FALSE]
}


# log2 and z-score normalization
norm_dat <- function(df, nor) {
  df <- log2(df+1)
  
  if (nor=="two.end") {
    df <- apply(df,1,function(x) (x-mean(x))/sd(x))
  }
  else if (nor=="up") 
    {
    df <- apply(df,1,function(x) ((x - min(x)) / (max(x) - min(x)))*(2-1)+1)
  }
  else if (nor=="down")
    {
    df <- apply(df, 1, function(x) ((max(x) - x) / (max(x) - min(x))) * (2 - 1) + 1)
  }
  else {print("wrong input for nor")}
  df <- t(df)
}


# making trees
make_tree <- function(dat, directory, group_label) {
  no_of_genes <- dim(dat)[1]
  gene_names <- rownames(dat)
  
  # Get the gene expression data
  ge_df_matrix <- as.matrix(dat)
  rownames(ge_df_matrix) <- rownames(dat)
  
  # Transpose the above matrix
  transposed_ge_df <- t(ge_df_matrix)
  transposed_ge_df <- as.data.frame(transposed_ge_df)
  
  
  # For each gene we make the tree
  for (gene_no in seq(no_of_genes)) {
    # Print the gene name and the number
    print(paste("Working on gene no.", gene_no, ":", gene_names[gene_no]))
    
    # Get the data for the gene
    gene_df <- as.data.frame(transposed_ge_df[gene_no])
    # head(gene_df)
    gene_name = colnames(gene_df)
    
    # Get the gene events
    ## Dissimilarity matrix based on euclidean distance
    d <- dist(gene_df, method = "euclidean")
    ## Hierarchical clustering using median linkage
    hc1 <- hclust(d, method = "median" )
    ## Generate dendrogram
    dend <- as.dendrogram(hc1)
    ## Get levels (depth)
    tree <- as.Node(dend)
    levels <- tree$Get("level")
    ## Return all the nodes in the tree with corresponding samples
    gene_events <- partition_leaves(dend)
    
    ## get dendrogram height
    height <- get_nodes_attr(dend, "height")
    
    
    # Create an empty data frame to append the result of each gene tree events
    df <- data.frame(matrix(ncol = 7, nrow = 0))
    
    for (i in seq(gene_events)) {
      
      # Get event characteristics
      event_no <- i
      # Get event samples
      unlisted_samples <- unlist(gene_events[i])
      samples <- paste(unlisted_samples, collapse=",")
      # Get median gene expression
      patients_ge <- gene_df %>% filter(row.names(gene_df) %in% unlisted_samples)
      median_exp <- median(patients_ge[[1]])
      
      # Get the number of samples and if there is only one sample, its a leaf
      no_of_samples <- length(unlisted_samples)
      if (no_of_samples == 1) {
        leaf_status <- 1
      } else {
        leaf_status <- 0
      }
      # Append the above data frame to df
      df <- rbind(df, c(event_no, samples, median_exp, no_of_samples, gene_name, 'gene_expression', leaf_status))
    }
    
    # add height
    df$hegith <- height
    df$levels <- levels
    # select top 10%
    df <- df[order(abs(as.numeric(df[,3])), decreasing = T),]
    nevent <- round(nrow(df) * 0.1)
    df$include <- c(rep("Y",nevent),rep("N",nrow(df)-nevent)) 
    df$group <- group_label
    # Provide column names
    colnames(df) <- c('event_no', 'samples', 'median_exp', 'no_of_samples', 'gene_name', 'event_name', 'leaf_status',"dendrogram_height",
                      "levels", "include", "group")
    
    file_name = paste(directory, gene_name, ".csv", sep="")
    write.table(df, file=file_name, sep=",",row.names = F)
  }
}


  return(p_fre_sub)
}

