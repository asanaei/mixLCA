# Compute Fit Indices for a mixLCA Model

Returns AIC, BIC, sample-adjusted BIC, relative entropy, and the
Integrated Classification Likelihood criterion (ICL).

## Usage

``` r
fit_indices(model)
```

## Arguments

- model:

  A `mixLCA` object.

## Value

Named list of fit statistics.

## Examples

``` r
# \donttest{
data(voter_perceptions)
fit <- fit_lca(voter_perceptions,
               categorical = names(voter_perceptions),
               n_classes   = 3,
               control     = lca_control(n_starts = 2),
               verbose     = FALSE)
fi <- fit_indices(fit)
fi[c("log_lik", "n_params", "AIC", "BIC", "aBIC", "entropy", "ICL")]
#> $log_lik
#> [1] -30210.13
#> 
#> $n_params
#> [1] 110
#> 
#> $AIC
#> [1] 60640.26
#> 
#> $BIC
#> [1] 61256.36
#> 
#> $aBIC
#> [1] 60906.89
#> 
#> $entropy
#> [1] 0.8210241
#> 
#> $ICL
#> [1] 62042.86
#> 
# }
```
