# Spectral Loadings as a Tidy Data Frame

Spectral Loadings as a Tidy Data Frame

## Usage

``` r
spectral_loadings(model)
```

## Arguments

- model:

  A `mixLCA` object.

## Value

Data frame.

## Examples

``` r
# \donttest{
data(voter_perceptions)
fit <- fit_lca(voter_perceptions,
               categorical   = names(voter_perceptions),
               n_classes     = 2,
               spectral_rank = c(1L, 1L),
               control       = lca_control(n_starts = 2, seed = 110),
               verbose       = FALSE)
ld <- spectral_loadings(fit)
head(ld, 12)
#>      class dimension     item  category     loading
#> 1  Class 1         1  moral_A excellent  0.27707145
#> 2  Class 1         1  moral_A      fair -0.14408882
#> 3  Class 1         1  moral_A      good -0.05305048
#> 4  Class 1         1  moral_A      poor -0.07993215
#> 5  Class 1         1  moral_B excellent  0.07210772
#> 6  Class 1         1  moral_B      fair  0.06958085
#> 7  Class 1         1  moral_B      good  0.07580412
#> 8  Class 1         1  moral_B      poor -0.21749269
#> 9  Class 1         1 compet_A excellent  0.27772845
#> 10 Class 1         1 compet_A      fair -0.09818432
#> 11 Class 1         1 compet_A      good -0.09227241
#> 12 Class 1         1 compet_A      poor -0.08727171
# }
```
