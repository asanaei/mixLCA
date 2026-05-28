# Print a mixDistal Object

Prints the per-class regression coefficients, sandwich standard errors,
and classification error matrix produced by
[`distal`](https://asanaei.github.io/mixLCA/reference/distal.md).

## Usage

``` r
# S3 method for class 'mixDistal'
print(x, ...)
```

## Arguments

- x:

  A `mixDistal` object.

- ...:

  Unused.

## Value

Invisibly returns `x`.

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
d <- distal(fit, health_screening, outcome ~ age, family = "binomial")
print(d)
#> 
#> Distal Outcome Estimation (BCH Method) - mixLCA
#> ================================================
#> Formula: outcome ~ age 
#> Family : binomial 
#> Classes: 2 
#> 
#> --- Class 1 ---
#>             Estimate       SE      z      p  
#> (Intercept) -1.67727  0.90877 -1.846 0.0649 .
#> age          0.02628  0.01766  1.488 0.1367  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>   Effective N: 167.3 
#> 
#> --- Class 2 ---
#>             Estimate      SE      z        p    
#> (Intercept)  -3.9106  0.7234 -5.406 6.44e-08 ***
#> age           0.0373  0.0154  2.422   0.0154 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>   Effective N: 632.7 
#> 
#> Classification Error Matrix:
#>        [,1]  [,2]
#> [1,] 0.8134 0.022
#> [2,] 0.1866 0.978
#> 
# }
```
