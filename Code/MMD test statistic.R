### compute MMD test statistic using kernel matrix


MMD_test_statistic <- function(kernel_matrix,n,m){
  
  diag(kernel_matrix) <- 0
  
  T_1 <- sum(kernel_matrix[1:n,1:n]) / (n^2 - n)
  
  T_2 <- sum(kernel_matrix[(n+1):(n+m),(n+1):(n+m)]) / (m^2 - m)
  
  T_3 <- 2*sum(kernel_matrix[1:n,(n+1):(n+m)])/(n*m)
  
  res <- T_1 + T_2 - T_3
  
  return(res)
}

