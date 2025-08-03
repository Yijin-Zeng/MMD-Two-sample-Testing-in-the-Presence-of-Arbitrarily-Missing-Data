# MMD two sample testing in the presence of arbitrarily missing data

The Maximum Mean Discrepancy (MMD) two-sample test (Gretton et al., 2012) is a popular non-parametric method capable of handling both univariate and multivariate data. When characteristic kernels—such as Gaussian or Laplacian—are employed, MMD can detect any distributional shift. However, it cannot be directly applied when the data contain missing values or are only partially observed.

This repository provides R implementations of Maximum Mean Discrepancy (MMD) two-sample tests designed to handle missing data, based on the methods proposed by (Zeng et al. 2024). These methods make no assumptions about the missingness mechanisms, and are particularly useful when data are missing not at random (MNAR). The methods are based on finding mathematical precise bounds of the MMD test statistics when the Laplacian kernel is used, and reject the null hypothesis when all possible test statistics are significant. For further details, see Zeng et al. (2024).

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

Gretton A, Borgwardt KM, Rasch MJ, Schölkopf B, Smola A. A kernel two-sample test. The journal of machine learning research. 2012 Mar 1;13(1):723-73.

Zeng Y, Adams NM, Bodenham DA. MMD Two-sample Testing in the Presence of Arbitrarily Missing Data. arXiv preprint arXiv:2405.15531. 2024 May 24.
