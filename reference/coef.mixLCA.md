# Coef Method for mixLCA

Returns a list of all estimated parameters.

## Usage

``` r
# S3 method for class 'mixLCA'
coef(object, ...)
```

## Arguments

- object:

  A `mixLCA` object.

- ...:

  Unused.

## Value

Named list.

## Examples

``` r
# \donttest{
data(health_screening)
fit <- fit_lca(health_screening,
               continuous  = c("marker_1","marker_2","marker_3","marker_4"),
               concomitant = ~ age,
               n_classes   = 2,
               control     = lca_control(n_starts = 2, seed = 110),
               verbose     = FALSE)
cc <- coef(fit)
names(cc)
#> [1] "concomitant" "means"       "covariances"
cc$concomitant
#>                 Class_2
#> (Intercept)  3.49083451
#> age         -0.04707317
cc$means
#> [[1]]
#> [1] 130.9492020  97.6234765  31.1183691   0.5408622
#> 
#> [[2]]
#> [1] 82.6076039 73.9968838 22.3790904  0.3120182
#> 
# }
```
