# MMD Testing with Missing Data

This repository contains R implementations for performing Maximum Mean Discrepancy (MMD) two-sample tests in the presence of missing data methods proposed in [1]. The methods provide both permutation-based and asymptotic testing procedures that can handle missing not at random (MNAR) data patterns.

## Overview

Maximum Mean Discrepancy (MMD) is a kernel-based statistical test used to determine whether two samples come from the same distribution. This implementation extends traditional MMD testing to handle scenarios where data contains missing values, which is common in real-world applications.

## Features

- **Two-sample testing with missing data**: Robust MMD testing when samples contain missing values
- **Multiple testing approaches**: 
  - Permutation-based testing (suitable for all sample sizes)
  - Asymptotic testing using Central Limit Theorem (recommended for n,m ≥ 25, d ≥ 50)
- **Univariate and multivariate support**: Handle both 1D and multi-dimensional data
- **Automatic bandwidth selection**: Median heuristic for kernel parameter selection
- **MCAR assumption**: Designed for missing completely at random data patterns

## Methods

The package implements bounds-based approaches:
- **Lower bounds** for MMD estimates using available data
- **Upper bounds** for variance estimation
- **Studentized test statistics** for improved finite-sample performance

## Files Structure

- `Examples.R` - Comprehensive examples demonstrating all testing scenarios
- `MMD using permutation with Missing data.R` - Permutation-based testing implementation
- `MMD using CLT with Missing data.R` - Asymptotic testing using CLT
- `Laplacian Kernel.R` - Kernel function implementations
- `Median Heuristic.R` - Automatic bandwidth selection
- `Univariate Bounds MMD.R` - Bounds computation for 1D data
- `Multivariate Lower Bounds MMD.R` - Lower bounds for multivariate data
- `Multivariate Upper Bounds MMD.R` - Upper bounds for multivariate data
- `Multivariate Maximum Variance MMD.R` - Variance estimation for multivariate data
- `CLT for studentized test statistic.R` - Asymptotic distribution theory

## Usage

### Basic Setup

```r
# Source the main functions
source('MMD using permutation with Missing data.R')
source('MMD using CLT with Missing data.R')

# Set random seed for reproducibility
set.seed(0)
```

### Example 1: Complete Univariate Data

```r
n <- 100
m <- 100
X <- rnorm(n, 0, 1)
Y <- rnorm(m, 0, 1)

# Automatic bandwidth selection
beta <- MedianHeuristic(X, Y)

# Perform permutation test
result <- permutation_testing_with_missing_data(X, Y, beta, perm = 200)
print(result$stat)  # Test statistic
print(result$pval)  # P-value
```

### Example 2: Univariate Data with Missing Values

```r
# Generate missing data (MCAR)
MCAR_univariate <- function(X, s) {
  n <- length(X)
  missing_location <- sample(1:n, s*n)
  X[missing_location] <- NA
  return(X)
}

n <- 100
m <- 100
X <- rnorm(n, 0, 1)
Y <- rnorm(m, 1, 1)  # Different means
s <- 0.05  # 5% missing data

# Introduce missing data
Incomplete_X <- MCAR_univariate(X, s)
Incomplete_Y <- MCAR_univariate(Y, s)

# Test with missing data
beta <- MedianHeuristic(Incomplete_X[!is.na(Incomplete_X)], 
                       Incomplete_Y[!is.na(Incomplete_Y)])
result <- permutation_testing_with_missing_data(Incomplete_X, Incomplete_Y, beta, perm = 100)
```

### Example 3: Multivariate Data with Missing Values

```r
library(MASS)  # Required for mvrnorm

# Generate multivariate data
d <- 10
n <- 100
m <- 100
mu_1 <- rep(0, d)
sigma_1 <- diag(d)
X <- mvrnorm(n, mu_1, sigma_1)
mu_2 <- rep(1, d)  # Different means
sigma_2 <- diag(d)
Y <- mvrnorm(m, mu_2, sigma_2)

# Introduce missing data pattern
MCAR_Multivariate <- function(X, S, s) {
  n <- dim(X)[1]
  d <- dim(X)[2]
  missing_location <- sample(1:n, S*n)
  X[missing_location, sample(1:d, s*d)] <- NA
  return(X)
}

S <- 0.05  # 5% of samples have missing values
s <- 0.2   # 20% of dimensions missing in incomplete samples
Incomplete_X <- MCAR_Multivariate(X, S, s)
Incomplete_Y <- MCAR_Multivariate(Y, S, s)

# Test with multivariate missing data
beta <- MedianHeuristic(Incomplete_X[!is.na(rowSums(Incomplete_X)),],
                       Incomplete_Y[!is.na(rowSums(Incomplete_Y)),])
result <- permutation_testing_with_missing_data(Incomplete_X, Incomplete_Y, beta, perm = 100)
```

### Example 4: Asymptotic Testing (High-Dimensional)

```r
# For high-dimensional data (d ≥ 50, n,m ≥ 25)
d <- 50
n <- 100
m <- 100
X <- mvrnorm(n, rep(0, d), diag(d))
Y <- mvrnorm(m, rep(0.8, d), diag(d))

# Introduce missing data
Incomplete_X <- MCAR_Multivariate(X, 0.05, 0.2)
Incomplete_Y <- MCAR_Multivariate(Y, 0.05, 0.2)

# Use CLT-based testing (faster for high dimensions)
beta <- MedianHeuristic(Incomplete_X[!is.na(rowSums(Incomplete_X)),],
                       Incomplete_Y[!is.na(rowSums(Incomplete_Y)),])
result <- testing_with_missing_using_CLT(Incomplete_X, Incomplete_Y, beta)
```

## Dependencies

- **Base R**: No external packages required for core functionality
- **MASS**: Required for multivariate normal random generation in examples
- **Optional**: `eummd` package for comparison with existing MMD implementations

## Method Selection Guidelines

- **Permutation testing**: Suitable for all sample sizes and dimensions. More computationally intensive but exact.
- **CLT-based testing**: Recommended when n,m ≥ 25 and d ≥ 50. Faster computation with asymptotic approximation.
- **Univariate**: Use permutation testing only, as normality approximation may not be suitable for low-dimensional data.

## Theory

This implementation is based on statistical theory for handling missing data in kernel-based two-sample testing. The methods compute:

1. **Lower bounds** on MMD estimates using available data
2. **Upper bounds** on variance estimates to maintain test validity
3. **Conservative p-values** that account for missing data uncertainty

## References

[1] Gretton A, Borgwardt KM, Rasch MJ, Schölkopf B, Smola A. A kernel two-sample test. The journal of machine learning research. 2012 Mar 1;13(1):723-73.

[2] Zeng Y, Adams NM, Bodenham DA. MMD Two-sample Testing in the Presence of Arbitrarily Missing Data. arXiv preprint arXiv:2405.15531. 2024 May 24.
