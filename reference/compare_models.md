# Compare Multiple mixLCA Models

Assembles a summary table of fit indices across a list of fitted models,
ordered by number of classes.

## Usage

``` r
compare_models(...)
```

## Arguments

- ...:

  One or more `mixLCA` objects, or a single named list of them.

## Value

Data frame of comparative fit statistics.

## Examples

``` r
# \donttest{
data(voter_perceptions)
cat_items <- names(voter_perceptions)
fits <- lapply(2:4, function(K)
  fit_lca(voter_perceptions, categorical = cat_items, n_classes = K,
          control = lca_control(n_starts = 2),
          verbose = FALSE))
names(fits) <- paste0("K", 2:4)
compare_models(fits)
#>    K        LL n_params      AIC      BIC     aBIC   entropy      ICL
#> K2 2 -30684.79       73 61515.59 61924.45 61692.53 0.8455514 62352.68
#> K3 3 -30210.13      110 60640.26 61256.36 60906.89 0.8210178 62042.89
#> K4 4 -29933.25      147 60160.50 60983.84 60516.81 0.7775981 62217.09
# }
```
