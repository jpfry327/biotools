# R/plotting.R

#' PCA plot of gene expression data
#' 
#' @param expr_matrix Gene expression matrix (genes as rows, samples as columns)
#' @param sample_labels Optional vector of sample labels
#' @param title Plot title
#' @param log_transform Whether to log2 transform the data
#' @return ggplot2 object
#' @export
#' @import ggplot2
#' @examples
#' \dontrun{
#' plot_gene_pca(expr_matrix)
#' }
plot_gene_pca <- function(expr_matrix, sample_labels = NULL, title = "PCA of Gene Expression", 
                         log_transform = FALSE) {
  # Remove genes with zero variance
  expr_clean <- expr_matrix[apply(expr_matrix, 1, var) > 0, ]
  
  # Optional log transformation
  if(log_transform) {
    expr_clean <- log2(expr_clean + 1)
  }
  
  # Transpose so samples are rows for PCA
  expr_t <- t(expr_clean)
  
  # Perform PCA
  pca_result <- prcomp(expr_t, scale. = TRUE, center = TRUE)
  
  # Extract variance explained
  var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
  
  # Create plot
  plot_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Sample = if(is.null(sample_labels)) rownames(pca_result$x) else sample_labels
  )
  
  ggplot(plot_data, aes(x = PC1, y = PC2, label = Sample)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, hjust = 0.5) +
    labs(
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)"),
      title = title
    ) +
    theme_minimal()
}

#' t-SNE plot of gene expression data
#' 
#' @param expr_matrix Gene expression matrix (genes as rows, samples as columns)
#' @param sample_labels Optional vector of sample labels
#' @param title Plot title
#' @param log_transform Whether to log2 transform the data
#' @param perplexity Perplexity parameter for t-SNE
#' @return ggplot2 object
#' @export
#' @import ggplot2
#' @import Rtsne
plot_gene_tsne <- function(expr_matrix, sample_labels = NULL, title = "t-SNE of Gene Expression", 
                          log_transform = FALSE, perplexity = 5) {
  # Remove genes with zero variance
  expr_clean <- expr_matrix[apply(expr_matrix, 1, var) > 0, ]
  
  # Optional log transformation
  if(log_transform) {
    expr_clean <- log2(expr_clean + 1)
  }
  
  # Transpose so samples are rows
  expr_t <- t(expr_clean)
  
  # Scale the data
  expr_scaled <- scale(expr_t)
  
  # Run t-SNE (adjust perplexity based on sample size)
  set.seed(42)  # For reproducibility
  tsne_result <- Rtsne(expr_scaled, dims = 2, perplexity = min(perplexity, (nrow(expr_scaled)-1)/3))
  
  # Create plot
  plot_data <- data.frame(
    tSNE1 = tsne_result$Y[, 1],
    tSNE2 = tsne_result$Y[, 2],
    Sample = if(is.null(sample_labels)) rownames(expr_t) else sample_labels
  )
  
  ggplot(plot_data, aes(x = tSNE1, y = tSNE2, label = Sample)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, hjust = 0.5) +
    labs(x = "t-SNE 1", y = "t-SNE 2", title = title) +
    theme_minimal()
}

#' MDS plot of gene expression data
#' 
#' @param expr_matrix Gene expression matrix (genes as rows, samples as columns)
#' @param sample_labels Optional vector of sample labels
#' @param title Plot title
#' @param log_transform Whether to log2 transform the data
#' @param distance_method Distance method for MDS
#' @return ggplot2 object
#' @export
#' @import ggplot2
plot_gene_mds <- function(expr_matrix, sample_labels = NULL, title = "MDS of Gene Expression", 
                         log_transform = FALSE, distance_method = "euclidean") {
  # Remove genes with zero variance
  expr_clean <- expr_matrix[apply(expr_matrix, 1, var) > 0, ]
  
  # Optional log transformation
  if(log_transform) {
    expr_clean <- log2(expr_clean + 1)
  }
  
  # Transpose so samples are rows
  expr_t <- t(expr_clean)
  
  # Scale the data
  expr_scaled <- scale(expr_t)
  
  # Calculate distance matrix
  dist_matrix <- dist(expr_scaled, method = distance_method)
  
  # Perform classical MDS
  mds_result <- cmdscale(dist_matrix, k = 2)
  
  # Create plot
  plot_data <- data.frame(
    MDS1 = mds_result[, 1],
    MDS2 = mds_result[, 2],
    Sample = if(is.null(sample_labels)) rownames(expr_t) else sample_labels
  )
  
  ggplot(plot_data, aes(x = MDS1, y = MDS2, label = Sample)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, hjust = 0.5) +
    labs(x = "MDS Dimension 1", y = "MDS Dimension 2", title = title) +
    theme_minimal()
}