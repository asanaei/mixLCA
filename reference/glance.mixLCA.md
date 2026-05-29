# Glance at a mixLCA Model

Returns a one-row data frame of fit statistics suitable for model
comparison tables and `tune_cluster()` logging.

## Usage

``` r
# S3 method for class 'mixLCA'
glance(x, ...)
```

## Arguments

- x:

  A `mixLCA` object.

- ...:

  Unused.

## Value

A one-row data frame with columns `n_classes`, `logLik`, `AIC`, `BIC`,
`aBIC`, `entropy`, `ICL`, `nobs`, `n_params`, `converged`, `iterations`.

## Examples

``` r
# \donttest{
data(voter_perceptions)
fit <- fit_lca(voter_perceptions,
               categorical = names(voter_perceptions),
               n_classes   = 3,
               control     = lca_control(n_starts = 2),
               verbose     = FALSE)
glance(fit)
#>   n_classes    logLik      AIC      BIC     aBIC   entropy      ICL nobs
#> 1         3 -30210.13 60640.26 61256.36 60906.89 0.8210178 62042.89 2000
#>   n_params converged iterations
#> 1      110      TRUE         68
# }
```
