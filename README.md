# MMD two sample testing in the presence of arbitrarily missing data

Maximum Mean Discrepancy (MMD) two-sample tests [1] are popular non-parametric testing methods that can handle both univariate and multivariate data. When characteristic kernels such as Gaussian or Laplacian kernels are used, MMD can detect any distributional shift. However, MMD can not be used directly when data may be missing, or only partially observed. 

This repository contains R implementations for performing Maximum Mean Discrepancy (MMD) two-sample tests in the presence of missing data methods proposed in [1]. The methods provide both permutation-based and asymptotic testing procedures that can handle missing not at random (MNAR) data patterns, which is common in real-world applications.

See more details in [1].

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

## References

[1] Gretton A, Borgwardt KM, Rasch MJ, Schölkopf B, Smola A. A kernel two-sample test. The journal of machine learning research. 2012 Mar 1;13(1):723-73.

[2] Zeng Y, Adams NM, Bodenham DA. MMD Two-sample Testing in the Presence of Arbitrarily Missing Data. arXiv preprint arXiv:2405.15531. 2024 May 24.
