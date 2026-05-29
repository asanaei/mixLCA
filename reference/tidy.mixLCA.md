# Tidy a mixLCA Model

Summarizes each latent class as one row, with size and mixing
proportion.

## Usage

``` r
# S3 method for class 'mixLCA'
tidy(x, ...)
```

## Arguments

- x:

  A `mixLCA` object.

- ...:

  Unused.

## Value

A data frame with columns `cluster`, `size`, `proportion`.

## Examples

``` r
# \donttest{
data(voter_perceptions)
fit <- fit_lca(voter_perceptions,
               categorical = names(voter_perceptions),
               n_classes   = 3,
               control     = lca_control(n_starts = 2),
               verbose     = FALSE)
tidy(fit)
#>     cluster size proportion
#> 1 Cluster_1  423  0.2114589
#> 2 Cluster_2  683  0.3413805
#> 3 Cluster_3  894  0.4471606
# }
```
