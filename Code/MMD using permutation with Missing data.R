### Compute bounds of MMD with missing data and performing test using permutations
source('Laplacian Kernel.R')
source('Median Heuristic.R')
source('Univariate Bounds MMD.R')
source('Multivariate Lower Bounds MMD.R')
source('Multivariate Upper Bounds MMD.R')


permutation_testing_with_missing_data <- function(X,Y,beta,perm){
  if(is.vector(X) & is.vector(Y)){
    n <- length(X)
    m <- length(Y)
    Z <- c(X,Y)
    test_statistic <- Lower_Bound_MMD_Univariate(X,Y,beta)
    res_upper_bound <- 0
    for (i in 1:perm) {
      permutation_X_Y <- sample(seq(1,(n+m)), (n+m))
      permutation_X <- Z[permutation_X_Y[1:(n)]]
      permutation_Y <- Z[permutation_X_Y[(n+1):(n+m)]]
      res_upper_bound[i] <- Upper_Bound_MMD_Univariate(permutation_X, permutation_Y, beta) 
    }
  }else{
    n <- dim(X)[1]
    m <- dim(Y)[1]
    Z<- rbind(X,Y)
    test_statistic <- Lower_Bound_MMD(X,Y,beta)
    res_upper_bound <- 0
    for (i in 1:perm) {
      permutation_X_Y <- sample(seq(1,(n+m)), (n+m) )
      permutation_X <- Z[permutation_X_Y[1:(n)],]
      permutation_Y <- Z[permutation_X_Y[(n+1):(n+m)],]
      res_upper_bound[i] <- Upper_Bound_MMD(permutation_X, permutation_Y, beta) 
    } 
  }
  p_value <- (sum(res_upper_bound >= test_statistic) +1)/(perm+1)
  return(list(stat = test_statistic, pval = p_value))
}