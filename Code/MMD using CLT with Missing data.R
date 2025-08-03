### Compute bounds of MMD and performing test
source('Laplacian Kernel.R')
source('Median Heuristic.R')
source('Multivariate Lower Bounds MMD.R')
source('CLT for studentized test statistic.R')
source('Multivariate Maximum Variance MMD.R')

testing_with_missing_using_CLT <- function(X,Y,beta){
  if(is.vector(X) | is.vector(Y)){
    return('Normality Approximation may not be suitable for univariate samples')
  }
  # compute the upper bounds of variance
  bounds_matrix <- compute_kernel_matrix_min_max(X,Y,beta)
  kernel_matrix_min <- bounds_matrix$kernel_matrix_min
  kernel_matrix_max <- bounds_matrix$kernel_matrix_max
  max_var <- var_max_estimator(kernel_matrix_min, kernel_matrix_max)
  # compute lower bounds of MMD
  n <- dim(X)[1]
  m <- dim(Y)[1]
  c_n_m <- 2/(n*(n-1)) + 4/(n*m) + 2/(m*(m-1))
  studendized_test_statistic <- Lower_Bound_MMD(X,Y,beta)/sqrt(c_n_m*max_var)
  # compute p_value
  if(studendized_test_statistic > 0){
    p_value <- pnorm(studendized_test_statistic,lower.tail = FALSE)
  }else{
    p_value <- 1
  }
  return(list(stat = studendized_test_statistic, pval = p_value))
}


# d <- 50
# n <- 50
# m <- 50
# sigma_1 <- diag(d)
# sigma_2 <- diag(d)
# X <- MASS::mvrnorm(n,rep(0,d),sigma_1)
# Y <- MASS::mvrnorm(n,c(rep(0,d*0.5),rep(1.5,d*0.5)),sigma_2)
# beta <- eummd::medianheuristic(rbind(X,Y))
# X[1:n*0.1,1:d*0.3] <- NA
# Y[1:m*0.1,1:d*0.3] <- NA
# testing_with_missing_using_CLT(X,Y,beta)
