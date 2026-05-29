# AIC Method for mixLCA

AIC Method for mixLCA

## Usage

``` r
# S3 method for class 'mixLCA'
AIC(object, ..., k = 2)
```

## Arguments

- object:

  A `mixLCA` object.

- ...:

  Unused.

- k:

  Numeric penalty per parameter (default 2).

## Value

Numeric AIC value.

## Examples

``` r
# \donttest{
data(voter_perceptions)
fit <- fit_lca(voter_perceptions,
               categorical = names(voter_perceptions),
               n_classes   = 2,
               control     = lca_control(n_starts = 2),
               verbose     = FALSE)
AIC(fit)
#> [1] 61515.59
# BIC-like penalty (k = log(N))
AIC(fit, k = log(nobs(fit)))
#> [1] 61924.45
# }
```
