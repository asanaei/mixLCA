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

## Examples

``` r
# \donttest{
data(health_screening)
fit <- fit_lca(health_screening,
               continuous = c("marker_1","marker_2","marker_3","marker_4"),
               n_classes  = 2,
               control    = lca_control(n_starts = 2),
               verbose    = FALSE)
bvr_tests(fit, health_screening)
#>       var1     var2 residual_cov       chi_sq   p_value
#> 1 marker_1 marker_2 0.1634600176 1.360841e-05 0.9970566
#> 2 marker_1 marker_3 0.0893743151 6.764780e-05 0.9934376
#> 3 marker_1 marker_4 0.0029215622 1.333686e-04 0.9907858
#> 4 marker_2 marker_3 0.0313518139 1.607196e-05 0.9968013
#> 5 marker_2 marker_4 0.0013511379 5.507291e-05 0.9940789
#> 6 marker_3 marker_4 0.0004971185 1.239660e-04 0.9911165
# }
```
