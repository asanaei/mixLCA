# logLik Method for mixLCA

logLik Method for mixLCA

## Usage

``` r
# S3 method for class 'mixLCA'
logLik(object, ...)
```

## Arguments

- object:

  A `mixLCA` object.

- ...:

  Unused.

## Value

Object of class `logLik`.

## Examples

``` r
# \donttest{
data(voter_perceptions)
fit <- fit_lca(voter_perceptions,
               categorical = names(voter_perceptions),
               n_classes   = 2,
               control     = lca_control(n_starts = 2, seed = 110),
               verbose     = FALSE)
logLik(fit)
#> 'log Lik.' -30684.79 (df=73)
# }
```
