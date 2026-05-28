# Standard Errors for Continuous Parameters (Numerical)

Computes approximate standard errors for class-specific means and
diagonal covariance elements via numerical second derivative of the
observed marginal log-likelihood. This respects Louis's Principle by
differentiating through the mixture, not the Q-function.

## Usage

``` r
continuous_se(model, data, step = 1e-04)
```

## Arguments

- model:

  A `mixLCA` object.

- data:

  Data frame used for estimation.

- step:

  Finite difference step size.

## Value

A list with elements `mean_se` (list of K named vectors) and `cov_se`
(list of K d x d matrices), or NULL if no continuous indicators were
specified.

## Examples

``` r
# \donttest{
data(health_screening)
fit <- fit_lca(health_screening,
               continuous = c("marker_1","marker_2","marker_3","marker_4"),
               n_classes  = 2,
               control    = lca_control(n_starts = 3, seed = 110),
               verbose    = FALSE)
ses <- continuous_se(fit, health_screening)
round(ses$mean_se[[1]], 4)   # SE of class-1 means
#> marker_1 marker_2 marker_3 marker_4 
#>   1.2427   1.0594   0.2444   0.0045 
# }
```
