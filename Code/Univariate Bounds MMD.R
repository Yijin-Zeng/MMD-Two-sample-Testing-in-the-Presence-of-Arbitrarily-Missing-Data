### Compute lower bounds of MMD for univariate samples (d = 1)
source('Laplacian Kernel.R')

T_1_Univariate <- function(X,beta){
  # pairwise distance in X
  return(sum(exp(-beta*matrix(dist(X, method = 'manhattan')))))
}

T_2_Univariate <- function(Y,beta){
  # pairwise distance in Y: should be the same as T_1
  return(sum(exp(-beta*matrix(dist(Y, method = 'manhattan')))))
}

T_3_Univariate <- function(X,Y,beta){
  # pairwise distance between X and Y: should be the same as T_1
  return(sum(exp(-beta*abs(outer(X, Y, "-")))))
}

compute_A_3_each_minimum_Univariate <- function(Observed_X, Observed_Y, beta, c_1, c_3){

  Observed_Z <- c(Observed_X, Observed_Y)
  positive_part <- c_1 * rowSums(exp(-beta*abs(outer(Observed_Z, Observed_X, '-'))))
  negative_part <- c_3 * rowSums(exp(-beta*abs(outer(Observed_Z, Observed_Y, '-'))))
  return(min(min(0,positive_part - negative_part)))
  
}

Lower_Bound_MMD_Univariate <- function(X,Y,beta){
  
  Observed_X <- X[!is.na(X)]
  Observed_Y <- Y[!is.na(Y)]
  n <- length(X)
  m <- length(Y)
  n_prime <- length(Observed_X)
  m_prime <- length(Observed_Y)
  
  c_1 <- 2/(n*(n-1))
  c_2 <- 2/(m*(m-1))
  c_3 <- 2/(n*m)

  ### compute termone
  A_1_1 <- 0
  A_1_2 <- 0
  A_1_3 <- c_3*(n-n_prime)*(m-m_prime)
  A_1 <- A_1_1 + A_1_2 - A_1_3
  
  ### compute termtwo
  A_2_1 <- c_1*T_1_Univariate(Observed_X,beta)
  A_2_2 <- c_2*T_2_Univariate(Observed_Y,beta)
  A_2_3 <- c_3*T_3_Univariate(Observed_X, Observed_Y, beta)
  A_2 <- A_2_1 + A_2_2 - A_2_3
  
  ### compute termthree
  A_3 <- 0
  if( (n - n_prime > 0) ){
    A_3 <- compute_A_3_each_minimum_Univariate(Observed_X, Observed_Y, beta, c_1, c_3)*(n-n_prime)
  }
  
  ### compute termfour
  A_4 <- 0
  if( (m - m_prime > 0) ){
    
    A_4 <- compute_A_3_each_minimum_Univariate(Observed_Y, Observed_X, beta, c_2, c_3)*(m-m_prime)
    
  }
  
  res <- A_1 + A_2 + A_3 + A_4
  return(res)
}


compute_A_3_each_maximum_Univariate <- function(Observed_X, Observed_Y, beta, c_1, c_3){
  
  Observed_Z <- c(Observed_X, Observed_Y)
  positive_part <- c_1 * rowSums(exp(-beta*abs(outer(Observed_Z, Observed_X, '-'))))
  negative_part <- c_3 * rowSums(exp(-beta*abs(outer(Observed_Z, Observed_Y, '-'))))
  return(max(max(0,positive_part - negative_part)))
  
}


Upper_Bound_MMD_Univariate <- function(X,Y,beta){
  
  Observed_X <- X[!is.na(X)]
  Observed_Y <- Y[!is.na(Y)]
  n <- length(X)
  m <- length(Y)
  n_prime <- length(Observed_X)
  m_prime <- length(Observed_Y)
  
  c_1 <- 2/(n*(n-1))
  c_2 <- 2/(m*(m-1))
  c_3 <- 2/(n*m)
  
  ### compute termone
  A_1_1 <- (n-n_prime)*(n-n_prime-1)*c_1/2
  A_1_2 <- (m-m_prime)*(m-m_prime-1)*c_2/2
  A_1_3 <- 0
  A_1 <- A_1_1 + A_1_2 - A_1_3
  
  ### compute termtwo
  A_2_1 <- c_1*T_1_Univariate(Observed_X,beta)
  A_2_2 <- c_2*T_2_Univariate(Observed_Y,beta)
  A_2_3 <- c_3*T_3_Univariate(Observed_X, Observed_Y, beta)
  A_2 <- A_2_1 + A_2_2 - A_2_3
  
  ### compute termthree
  A_3 <- 0
  if( (n - n_prime > 0) ){
    A_3 <- compute_A_3_each_maximum_Univariate(Observed_X, Observed_Y, beta, c_1, c_3)*(n-n_prime)
  }
  
  ### compute termfour
  A_4 <- 0
  if( (m - m_prime > 0) ){
    A_4 <- compute_A_3_each_maximum_Univariate(Observed_Y, Observed_X, beta, c_2, c_3)*(m-m_prime)
  }
  
  res <- A_1 + A_2 + A_3 + A_4
  return(res)
}




