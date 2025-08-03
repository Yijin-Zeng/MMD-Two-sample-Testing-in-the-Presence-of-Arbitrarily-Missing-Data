############## Compute $p$-value of proposed method in the presence missing data.
rm(list = ls())
set.seed(0)
source('MMD using permutation with Missing data.R')
source('MMD using CLT with Missing data.R')

####### For univariate samples (Proposed: Perm method only, since normality 
####### approximation may only be suitable for n,m >=25, d >= 50).

## case 1: samples are all observed
n <- 100
m <- 100
X <- rnorm(n, 0, 1)
Y <- rnorm(m, 0, 1)
beta <- MedianHeuristic(X,Y)
res_case_1 <- permutation_testing_with_missing_data(X,Y,beta,perm = 200)

res_case_1$stat
res_case_1$pval

#### (optional) compare with results from eummd package
# res_eummd <- eummd::eummd(X,Y)
# res_eummd$stat  
# res_eummd$pval

## case 2: samples with missing data
MCAR_univariate <- function(X,s){
  # given X, return incomplete X with s proportion of missing data
  n <- length(X)
  missing_location <- sample(1:n, s*n)
  X[missing_location] <- NA
  return(X)
}

n <- 100
m <- 100
X <- rnorm(n, 0, 1)
Y <- rnorm(m, 1, 1)
s <- 0.05
Incomplete_X <- MCAR_univariate(X, s)
Incomplete_Y <- MCAR_univariate(Y, s)
beta <- MedianHeuristic(Incomplete_X[!is.na(Incomplete_X)], Incomplete_Y[!is.na(Incomplete_Y)])
res_case_2 <- permutation_testing_with_missing_data(Incomplete_X,Incomplete_Y,beta,perm = 100)
res_case_2$stat
res_case_2$pval  

##### For multivariate samples based on permutation (Proposed: Perm)
##### MASS packages are needed for generating random samples 
##### from multi-normal distribution

## case 3:  samples are all observed
d <- 10
n <- 100
m <- 100
mu_1 <- rep(0,d)
sigma_1 <- diag(d)
X <- MASS::mvrnorm(n, mu_1, sigma_1)
mu_2 <- rep(0,d)
sigma_2 <- diag(d)
Y <- MASS::mvrnorm(n, mu_2, sigma_2)
beta <- MedianHeuristic(X,Y)
res_case_3 <- permutation_testing_with_missing_data(X, Y, beta, perm = 100)
res_case_3$stat
res_case_3$pval

#### (optional) compare with results from eummd package
#res_case_3_eummd <- eummd::mmd(X,Y)
#res_case_3_eummd$stat  
#res_case_3_eummd$pval

## case 4:  not all samples are completely observed
MCAR_Multivariate <- function(X,S,s){
  # given X, return incomplete X with S proportion of incompletely observed samples,
  # each incomplete sample with s proportion of unobserved dimensions 
  n <- dim(X)[1]
  d <- dim(X)[2]
  missing_location <- sample(1:n, S*n)
  X[missing_location, sample(1:d, s*d)] <- NA
  return(X)
}

d <- 10
n <- 100
m <- 100
mu_1 <- rep(0,d)
sigma_1 <- diag(d)
X <- MASS::mvrnorm(n, mu_1, sigma_1)
mu_2 <- rep(1,d)
sigma_2 <- diag(d)
Y <- MASS::mvrnorm(n, mu_2, sigma_2)
S <- 0.05
s <- 0.2
Incomplete_X <- MCAR_Multivariate(X,S,s)
Incomplete_Y <- MCAR_Multivariate(Y,S,s)

beta <- MedianHeuristic(Incomplete_X[!is.na(rowSums(Incomplete_X)),],Incomplete_Y[!is.na(rowSums(Incomplete_Y)),])
### This may takes a few seconds
res_case_4 <- permutation_testing_with_missing_data(Incomplete_X, Incomplete_Y, beta, perm = 100)
res_case_4$stat
res_case_4$pval

##### For multivariate samples based on normality approximation (Proposed: Normality), 
##### when n,m >= 25, d >= 50


## case 5:  samples are all observed
d <- 50
n <- 100
m <- 100
mu_1 <- rep(0,d)
sigma_1 <- diag(d)
X <- MASS::mvrnorm(n, mu_1, sigma_1)
mu_2 <- rep(0,d)
sigma_2 <- diag(d)
Y <- MASS::mvrnorm(n, mu_2, sigma_2)
beta <- MedianHeuristic(X,Y)
res_case_5 <- testing_with_missing_using_CLT(X, Y, beta)
res_case_5$stat
res_case_5$pval

## case 6:  not all samples are completely observed

d <- 50
n <- 100
m <- 100
mu_1 <- rep(0,d)
sigma_1 <- diag(d)
X <- MASS::mvrnorm(n, mu_1, sigma_1)
mu_2 <- rep(0.8,d)
sigma_2 <- diag(d)
Y <- MASS::mvrnorm(n, mu_2, sigma_2)
S <- 0.05
s <- 0.2
Incomplete_X <- MCAR_Multivariate(X,S,s)
Incomplete_Y <- MCAR_Multivariate(Y,S,s)

beta <- MedianHeuristic(Incomplete_X[!is.na(rowSums(Incomplete_X)),],Incomplete_Y[!is.na(rowSums(Incomplete_Y)),])
res_case_6 <- testing_with_missing_using_CLT(Incomplete_X, Incomplete_Y, beta)
res_case_6$stat
res_case_6$pval

