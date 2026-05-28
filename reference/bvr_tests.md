# Bivariate Residual Significance Tests

Converts residual covariances to approximate chi-squared statistics (df
= 1) for pairwise local dependence testing.

## Usage

``` r
bvr_tests(model, data)
```

## Arguments

- model:

  A `mixLCA` object.

- data:

  Data frame used for estimation.

## Value

Data frame with columns: `var1`, `var2`, `residual_cov`, `chi_sq`,
`p_value`.
