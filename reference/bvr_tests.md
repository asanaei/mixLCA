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
#> 1 marker_1 marker_2 0.1634123576 1.360048e-05 0.9970575
#> 2 marker_1 marker_3 0.0893441432 6.760213e-05 0.9934398
#> 3 marker_1 marker_4 0.0029205507 1.332763e-04 0.9907890
#> 4 marker_2 marker_3 0.0313433284 1.606326e-05 0.9968022
#> 5 marker_2 marker_4 0.0013506172 5.503047e-05 0.9940811
#> 6 marker_3 marker_4 0.0004969341 1.238740e-04 0.9911198
# }
```
