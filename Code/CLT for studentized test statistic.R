source('Laplacian Kernel.R')
source('Median Heuristic.R')
source('MMD test statistic.R')

compute_kernel_matrix <- function(X, Y, beta) {
  
  Z <- rbind(X,Y)
  return(exp(-beta*as.matrix(dist(Z, method = "manhattan"))))
}


var_estimator <- function(kernel_matrix) {
  N <- dim(kernel_matrix)[1]
  
  # Compute a_dot_dot
  a_dot_dot <- sum(kernel_matrix) / ((N-1) * (N-2))
  
  # Compute a_dot_t and a_s_dot
  row_sums <- rowSums(kernel_matrix)
  col_sums <- colSums(kernel_matrix)
  a_dot_t <- col_sums / (N-2)
  a_s_dot <- row_sums / (N-2)
  
  # Compute A_s_t
  A_s_t <- kernel_matrix - matrix(rep(a_dot_t, each = N), ncol = N, byrow = FALSE) - matrix(rep(a_s_dot, each = N), nrow = N, byrow = TRUE) + a_dot_dot
  
  # Zero out the diagonal
  diag(A_s_t) <- 0
  
  # Compute the result
  result <- sum(A_s_t^2) / (N * (N-3)) - 1 / ((N-1) * (N-3))
  
  return(result)
}


testing_using_CLT <- function(X,Y,beta){
  ### computing p-value of MMD test statistic using normality approximation
  if( is.vector(X) | is.vector(Y) ){
    return('Normality Approximation may not be suitable for univariate samples')
  }
  n <- dim(X)[1]
  m <- dim(Y)[1]
  d <- dim(X)[2]
  c_n_m <- 2/(n*(n-1)) + 4/(n*m) + 2/(m*(m-1))
  kernel_matrix <- compute_kernel_matrix(X,Y,beta)
  estimated_variance <- var_estimator(kernel_matrix)
  test_statistic <- MMD_test_statistic(kernel_matrix,n,m)
  studenized_test_statistic <- test_statistic/sqrt(c_n_m*estimated_variance)
  p_value <- pnorm(studenized_test_statistic,lower.tail = FALSE)
  return(list(stat = studenized_test_statistic, pval = p_value))
}

# set.seed(1)
# d <- 50
# sigma_1 <- diag(d)
# sigma_2 <- diag(d)
# num <- 1e3
# dis_null <- rep(0,num)
# dis_pval <- rep(0,num)
# for(i in 1:num){
#   print(i)
#   X <- MASS::mvrnorm(500,rep(0,d),sigma_1)
#   Y <- MASS::mvrnorm(500,rep(0,d),sigma_2)
#   beta <- MedianHeuristic(X,Y)
#   testing_res <- testing_using_CLT(X,Y,beta) 
#   dis_null[i] <- testing_res$stat
#   dis_p_val[i] <- testing_res$pval
# }
# hist(dis_null)
# hist(dis_p_val)
