# Classification Table for a mixLCA Model

Cross-tabulates modal class assignments with average posterior
probabilities. Rows = modal assignment, columns = average P(class k).
Each row sums to 1. Diagonal dominance indicates good class separation.

## Usage

``` r
class_table(model)
```

## Arguments

- model:

  A `mixLCA` object.

## Value

K x K matrix: rows = modal assignment, columns = average P(class k).

## Examples

``` r
# \donttest{
data(voter_perceptions)
fit <- fit_lca(voter_perceptions,
               categorical = names(voter_perceptions),
               n_classes   = 3,
               control     = lca_control(n_starts = 2, seed = 110),
               verbose     = FALSE)
round(class_table(fit), 3)
#>       [,1]  [,2]  [,3]
#> [1,] 0.859 0.056 0.085
#> [2,] 0.024 0.957 0.019
#> [3,] 0.058 0.028 0.915
# }
```
