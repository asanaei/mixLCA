# Extract Posterior Class Probabilities

Accessor returning the N x K posterior probability matrix. Use this
rather than `model$posteriors` so downstream code is insulated against
future restructuring of internal fields.

## Usage

``` r
get_posteriors(x, ...)

# S3 method for class 'mixLCA'
get_posteriors(x, ...)
```

## Arguments

- x:

  A fitted `mixLCA` object.

- ...:

  Unused.

## Value

Numeric matrix.

## Examples

``` r
# \donttest{
data(voter_perceptions)
fit <- fit_lca(voter_perceptions,
               categorical = names(voter_perceptions),
               n_classes   = 2,
               control     = lca_control(n_starts = 2),
               verbose     = FALSE)
post <- get_posteriors(fit)
dim(post)
#> [1] 2000    2
head(round(post, 3))
#>       [,1]  [,2]
#> [1,] 0.983 0.017
#> [2,] 0.000 1.000
#> [3,] 0.997 0.003
#> [4,] 1.000 0.000
#> [5,] 1.000 0.000
#> [6,] 0.347 0.653
colSums(post) / nrow(post)   # estimated class proportions
#> [1] 0.4847926 0.5152074
# }
```
