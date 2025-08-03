source('Laplacian Kernel.R')
source('Multivariate Lower Bounds MMD.R')

IncompleteLaplacianKernel <- function(X,Y,beta){
  Z <- X - Y 
  return(exp(-beta*sum(abs(Z[!is.na(Z)]))))
}

compute_kernel_matrix_min_max <- function(X, Y, beta) {
  n <- nrow(X)
  m <- nrow(Y)  
  N <- n + m 
  Z <- rbind(X, Y)
  kernel_matrix_max <- matrix(NA, nrow = N, ncol = N)
  kernel_matrix_min <- matrix(NA, nrow = N, ncol = N)  
  
  # Samples completely observed
  complete_samples_X_indices <- which(!rowSums(is.na(X)))
  complete_samples_Y_indices <- which(!rowSums(is.na(Y)))
  
  complete_samples <- !rowSums(is.na(Z))
  
  # Compute pairwise distances for completely observed samples
  complete_indices <- which(complete_samples)
  complete_Z <- Z[complete_indices,]
  
  dist_observed <- exp(-beta * as.matrix(dist(complete_Z, method = "manhattan")))
  
  # Fill the kernel matrices for completely observed samples
  kernel_matrix_max[complete_indices, complete_indices] <- dist_observed
  kernel_matrix_min[complete_indices, complete_indices] <- dist_observed
  
  # Handle pairs with at least one sample having missing values
  for (i in 1:N) {
    for (j in 1:i) {
      if (is.na(kernel_matrix_max[i, j])) {
          dist_value <- IncompleteLaplacianKernel(Z[i, ], Z[j, ], beta)
          kernel_matrix_max[i, j] <- dist_value
          kernel_matrix_max[j, i] <- dist_value
          kernel_matrix_min[i, j] <- 0
          kernel_matrix_min[j, i] <- 0
      }
    }
  }
  diag(kernel_matrix_min) <- 1
  return(list(kernel_matrix_max = kernel_matrix_max, kernel_matrix_min = kernel_matrix_min))

}


var_max_estimator <- function(kernel_matrix_min, kernel_matrix_max){
  N <- dim(kernel_matrix_min)[1]
  
  # Compute a_dot_dot
  a_dot_dot <- sum(kernel_matrix_max) / ((N-1) * (N-2))
  
  # Compute a_dot_t and a_s_dot
  row_sums <- rowSums(kernel_matrix_min)
  col_sums <- colSums(kernel_matrix_min)
  a_dot_t <- col_sums / (N-2)
  a_s_dot <- row_sums / (N-2)
  
  # Compute A_s_t
  A_s_t <- kernel_matrix_max - matrix(rep(a_dot_t, each = N), ncol = N, byrow = FALSE) - matrix(rep(a_s_dot, each = N), nrow = N, byrow = TRUE) + a_dot_dot
  
  # Zero out the diagonal
  diag(A_s_t) <- 0
  
  # Compute the result
  result <- sum(A_s_t^2) / (N * (N-3)) - 1 / ((N-1) * (N-3))
  
  return(result)
}



