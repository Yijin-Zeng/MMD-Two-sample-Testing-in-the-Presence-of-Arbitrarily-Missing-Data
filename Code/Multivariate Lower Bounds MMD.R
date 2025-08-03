### Compute lower bounds of MMD for multivariate samples (d > 1)
source('Laplacian Kernel.R')
source('Median Heuristic.R')

T_1 <- function(X,beta){
  # pairwise distance in X
  return(sum(exp(-beta*matrix(dist(X, method = 'manhattan')))))
}

T_2 <- function(Y,beta){
  # pairwise distance in Y: should be the same as T_1
  return(sum(exp(-beta*matrix(dist(Y, method = 'manhattan')))))
}

T_3 <- function(X,Y,beta){
  # pairwise distance between X and Y: should be the same as T_1
  n <- dim(X)[1]
  m <- dim(Y)[1]
  d <- dim(X)[2]
  Res_T_3 <- 0
  ### compute T_3
  D_XY <- matrix(0, n, m)
  for (k in 1:d) {
    D_XY <- D_XY + abs(outer(X[, k], Y[, k], "-"))
  }
  return(sum(exp(-beta*D_XY)))
}

compute_A_3_each_minimum <- function(Incomplete_Vector_X, Complete_X, Complete_Y, beta, c_1, c_3){
  ### return the minimum value of each function of A_3
  
  d <- length(Incomplete_Vector_X)
  n_prime <- dim(Complete_X)[1]
  m_prime <- dim(Complete_Y)[1]
  missing_dimensions <- which(is.na(Incomplete_Vector_X))
  observed_dimensions <- setdiff(1:d, missing_dimensions)
  missing_d <- length(missing_dimensions)
  
  ### compute positive part
  Complete_Z <- rbind(Complete_X, Complete_Y)
  ### compute tilde x_i_j
  tilde_x_i_j <- matrix(0, nrow = n_prime, ncol = missing_d)
  for (i in seq(1,n_prime)) {
    if(length(missing_dimensions) > 1){
      max_dis_samples <- max.col(t(abs(sweep(Complete_Z[,missing_dimensions],2, Complete_X[i,missing_dimensions]))))
      tilde_x_i_j[i,] <- abs(Complete_X[i,missing_dimensions] - diag(Complete_Z[max_dis_samples,missing_dimensions]) ) 
    }else{
      tilde_x_i_j[i,1] <- max(abs(Complete_Z[,missing_dimensions] - Complete_X[i,missing_dimensions]))
    }
    
    
  }
  
  D_XX <- matrix(0, 1, n_prime)
  
  for (k in observed_dimensions ) {
    D_XX <- D_XX + abs(outer(Incomplete_Vector_X[k], Complete_X[,k], "-"))
  }
  
  positive_part <- c_1*sum(exp(-beta*D_XX)*exp(-beta*rowSums(tilde_x_i_j)))
  
  ### compute nagative part
  D_XY <- matrix(0, 1, m_prime)
  
  for (k in observed_dimensions) {
    D_XY <- D_XY + abs(outer(Incomplete_Vector_X[k], Complete_Y[,k], "-"))
  }
  
  negative_part <- c_3*sum(exp(-beta*D_XY))
  
  return(min(0, positive_part - negative_part))
}



compute_A_1_3 <- function(Incomplete_X, Incomplete_Y,beta){
  
  if(is.vector(Incomplete_X)){
    Incomplete_X <- matrix(Incomplete_X, nrow = 1)
  }
  
  if(is.vector(Incomplete_Y)){
    Incomplete_Y <- matrix(Incomplete_Y, nrow = 1)
  }
  
  n <- dim(Incomplete_X)[1]
  m <- dim(Incomplete_Y)[1]
  d <- dim(Incomplete_X)[2]
  
  D_XY <- matrix(0, n, m)
  
  for (k in 1:d) {
    D_XY_k <- abs(outer(Incomplete_X[, k], Incomplete_Y[, k], "-"))
    D_XY_k[is.na(D_XY_k)] <- 0
    D_XY <- D_XY + D_XY_k
  }
  
  return(sum(exp(-beta*D_XY)))

}

Lower_Bound_MMD <- function(X,Y,beta){
  
  n <- dim(X)[1]
  m <- dim(Y)[1]
  Complete_X <- X[! rowSums(is.na(X)) > 0,] ### samples in X completely observed
  Complete_Y <- Y[! rowSums(is.na(Y)) > 0,]  
  Incomplete_X <- X[rowSums(is.na(X)) > 0,] ### samples in X incompletely observed
  Incomplete_Y <- Y[rowSums(is.na(Y)) > 0,]
  n_prime <- dim(Complete_X)[1]
  m_prime <- dim(Complete_Y)[1]
  
  c_1 <- 2/(n*(n-1))
  c_2 <- 2/(m*(m-1))
  c_3 <- 2/(n*m)
  
  ### compute termone
  A_1_1 <- 0
  A_1_2 <- 0
  if( (n == n_prime) | (m == m_prime) ){
    A_1_3 <- 0
  }else{
    A_1_3 <- c_3*compute_A_1_3(Incomplete_X, Incomplete_Y,beta)
  }
  A_1 <- A_1_1 + A_1_2 - A_1_3
  
  ### compute termtwo
  A_2_1 <- c_1*T_1(Complete_X,beta)
  A_2_2 <- c_2*T_2(Complete_Y,beta)
  A_2_3 <- c_3*T_3(Complete_X, Complete_Y, beta)
  A_2 <- A_2_1 + A_2_2 - A_2_3
  
  
  ### compute termthree
  A_3 <- 0
  if( (n - n_prime > 0) ){
    
    
    if((n - n_prime == 1)){
      Incomplete_X <- matrix(Incomplete_X, nrow = 1)
    }
    
    for (i in 1:(n-n_prime)) {
      A_3 <- A_3 + compute_A_3_each_minimum(Incomplete_X[i,],Complete_X,Complete_Y,beta, c_1, c_3)
    }
    
  }
    

  ### compute termfour
  A_4 <- 0
  if( (m - m_prime > 0) ){
    
    
    if((m - m_prime == 1)){
      Incomplete_Y <- matrix(Incomplete_Y, nrow = 1)
    }
    
    for (i in 1:(m-m_prime)) {
      A_4 <- A_4 + compute_A_3_each_minimum(Incomplete_Y[i,],Complete_Y,Complete_X,beta,c_2, c_3)
    }
    
  }
  
  res <- A_1 + A_2 + A_3 + A_4
  return(res)
}

