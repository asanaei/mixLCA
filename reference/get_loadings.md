# Extract Spectral Loadings

Accessor returning the SLD loadings table (or NULL when the model does
not use Spectral Local Dependence). Use this rather than
`model$cat_spectral_params` so downstream code is insulated against
future restructuring of internal fields.

## Usage

``` r
get_loadings(x, ...)

# S3 method for class 'mixLCA'
get_loadings(x, ...)
```

## Arguments

- x:

  A fitted `mixLCA` object.

- ...:

  Unused.

## Value

A data frame or NULL.

## Examples

``` r
# \donttest{
data(voter_perceptions)
fit <- fit_lca(voter_perceptions,
               categorical   = names(voter_perceptions),
               n_classes     = 2,
               spectral_rank = c(1L, 1L),
               control       = lca_control(n_starts = 2),
               verbose       = FALSE)
ld <- get_loadings(fit)
head(ld, 12)
#>      class dimension     item  category     loading
#> 1  Class 1         1  moral_A excellent -0.05736222
#> 2  Class 1         1  moral_A      fair -0.08992328
#> 3  Class 1         1  moral_A      good -0.07013618
#> 4  Class 1         1  moral_A      poor  0.21742168
#> 5  Class 1         1  moral_B excellent -0.27767245
#> 6  Class 1         1  moral_B      fair  0.12489919
#> 7  Class 1         1  moral_B      good  0.04322520
#> 8  Class 1         1  moral_B      poor  0.10954806
#> 9  Class 1         1 compet_A excellent -0.05776005
#> 10 Class 1         1 compet_A      fair -0.07875409
#> 11 Class 1         1 compet_A      good -0.06689978
#> 12 Class 1         1 compet_A      poor  0.20341393
# }
```
