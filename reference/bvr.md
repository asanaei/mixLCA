# Bivariate Residual Covariance Matrix

Computes the residual covariance between continuous indicators after
subtracting the model-implied mixture covariance. Large residuals signal
local dependence violations not captured by the current model. Use
`plot_bvr()` to visualise the residual network.

## Usage

``` r
bvr(model, data)
```

## Arguments

- model:

  A `mixLCA` object.

- data:

  Data frame used for estimation.

## Value

Named matrix of residual covariances among continuous indicators.

## Examples

``` r
# \donttest{
data(health_screening)
fit <- fit_lca(health_screening,
               continuous = c("marker_1","marker_2","marker_3","marker_4"),
               n_classes  = 2,
               control    = lca_control(n_starts = 2, seed = 110),
               verbose    = FALSE)
round(bvr(fit, health_screening), 3)
#>          marker_1 marker_2 marker_3 marker_4
#> marker_1    2.197    0.163    0.089    0.003
#> marker_2    0.163    1.135    0.031    0.001
#> marker_3    0.089    0.031    0.068    0.000
#> marker_4    0.003    0.001    0.000    0.000
# }
```
