# Evaluate Continuous Log-Density

Computes the log-density of each observation under a class-specific
multivariate normal distribution. When entries are missing the density
is evaluated over the marginal distribution of the observed subset,
which constitutes exact analytical marginalisation rather than listwise
deletion or ad hoc imputation.

## Usage

``` r
eval_continuous_density(Y, mu, Sigma)
```

## Arguments

- Y:

  Numeric matrix (N x d), may contain `NA`.

- mu:

  Numeric vector of length d (class-specific mean).

- Sigma:

  Numeric d x d covariance matrix.

## Value

Numeric vector of length N (log-densities).
