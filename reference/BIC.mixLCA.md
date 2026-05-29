# BIC Method for mixLCA

BIC Method for mixLCA

## Usage

``` r
# S3 method for class 'mixLCA'
BIC(object, ...)
```

## Arguments

- object:

  A `mixLCA` object.

- ...:

  Unused.

## Value

Numeric BIC value.

## Examples

``` r
# \donttest{
data(voter_perceptions)
fit <- fit_lca(voter_perceptions,
               categorical = names(voter_perceptions),
               n_classes   = 2,
               control     = lca_control(n_starts = 2),
               verbose     = FALSE)
BIC(fit)
#> [1] 61924.45
# }
```
