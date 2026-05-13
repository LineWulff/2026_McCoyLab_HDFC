find_high_cor_pairs <- function(cor_mat, threshold = 0.5) {
  
  # Basic checks
  if (!is.matrix(cor_mat)) stop("Input must be a matrix")
  if (is.null(rownames(cor_mat)) || is.null(colnames(cor_mat))) {
    stop("Matrix must have row and column names")
  }
  
  # Get upper triangle indices (avoid duplicates + diagonal)
  upper_idx <- upper.tri(cor_mat, diag = FALSE)
  
  # Extract values
  vals <- cor_mat[upper_idx]
  
  # Get row/col names for those positions
  row_names <- rownames(cor_mat)[row(cor_mat)[upper_idx]]
  col_names <- colnames(cor_mat)[col(cor_mat)[upper_idx]]
  
  # Filter by threshold (exclude 1s automatically since diag removed)
  keep <- vals > threshold
  
  # Return as data frame
  data.frame(
    feature1 = row_names[keep],
    feature2 = col_names[keep],
    correlation = vals[keep],
    row.names = NULL
  )
}