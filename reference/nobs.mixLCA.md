# nobs Method for mixLCA

nobs Method for mixLCA

## Usage

``` r
# S3 method for class 'mixLCA'
nobs(object, ...)
```

## Arguments

- object:

  A `mixLCA` object.

- ...:

  Unused.

## Value

Integer.

## Examples

``` r
# \donttest{
data(voter_perceptions)
fit <- fit_lca(voter_perceptions,
               categorical = names(voter_perceptions),
               n_classes   = 2,
               control     = lca_control(n_starts = 2),
               verbose     = FALSE)
nobs(fit)
#> [1] 2000
# }
```
