# Print a mixLCA Object

Print a mixLCA Object

## Usage

``` r
# S3 method for class 'mixLCA'
print(x, ...)
```

## Arguments

- x:

  A `mixLCA` object.

- ...:

  Unused.

## Value

Invisibly returns `x`.

## Examples

``` r
# \donttest{
data(voter_perceptions)
fit <- fit_lca(voter_perceptions,
               categorical = names(voter_perceptions),
               n_classes   = 2,
               control     = lca_control(n_starts = 2),
               verbose     = FALSE)
print(fit)
#> 
#> Latent Class Model - mixLCA
#> ===========================
#> Classes        : 2 
#> Log-likelihood       : -30684.79 
#> Parameters     : 73 
#> Observations   : 2000 
#> Converged      : TRUE (27 iterations) 
#> 
#> Specification:
#>   Categorical indicators: 12 
#>   Dependence structure  : full 
#> 
#> Class proportions:
#> Class 1 Class 2 
#>  0.4848  0.5152 
#> 
# }
```
