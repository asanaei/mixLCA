# Fit Bridge for tidyclust

Translates the tidyclust `fit()` call into a
[`fit_lca()`](https://asanaei.github.io/mixLCA/reference/fit_lca.md)
call. Column roles (continuous vs. categorical) are specified through
engine arguments; when omitted they are inferred from column types.
Concomitant variables are excluded from auto-detection to prevent
double-counting.

## Usage

``` r
.lca_clust_fit_mixLCA(
  x,
  num_clusters,
  dependence = "full",
  penalty = NULL,
  spectral_rank = NULL,
  ...
)
```

## Arguments

- x:

  Data frame of predictors (intercept already removed by tidyclust
  encoding).

- num_clusters:

  Integer number of latent classes.

- dependence:

  Covariance structure passed to
  [`fit_lca()`](https://asanaei.github.io/mixLCA/reference/fit_lca.md).

- ...:

  Additional engine arguments forwarded to
  [`fit_lca()`](https://asanaei.github.io/mixLCA/reference/fit_lca.md)
  (e.g. `categorical`, `continuous`, `concomitant`, `control`).

## Value

A fitted `mixLCA` object.
