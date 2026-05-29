# Summary Method for mixLCA

Prints fit indices, class proportions, measurement parameters, and
concomitant coefficients with standard errors when `data` is supplied.

## Usage

``` r
# S3 method for class 'mixLCA'
summary(object, data = NULL, ...)
```

## Arguments

- object:

  A `mixLCA` object.

- data:

  Optional data frame: required for concomitant standard errors.

- ...:

  Unused.

## Value

Invisibly returns `object`.

## Examples

``` r
# \donttest{
data(health_screening)
fit <- fit_lca(health_screening,
               continuous  = c("marker_1","marker_2","marker_3","marker_4"),
               concomitant = ~ age,
               n_classes   = 2,
               control     = lca_control(n_starts = 3),
               verbose     = FALSE)
summary(fit, data = health_screening)
#> 
#> Summary - mixLCA
#> ================
#> 
#> Fit Indices:
#>   Log-likelihood     : -9938.795 
#>   AIC               : 19937.59 
#>   BIC               : 20078.13 
#>   Sample-adj. BIC   : 19982.86 
#>   Entropy           : 0.783 
#>   ICL               : 20318.83 
#> 
#> Class Proportions:
#> Class 1 Class 2 
#>  0.2092  0.7908 
#> 
#> Continuous Indicator Means:
#>         marker_1 marker_2 marker_3 marker_4
#> Class 1 130.9492  97.6235  31.1184   0.5409
#> Class 2  82.6076  73.9969  22.3791   0.3120
#> 
#> Covariance (Class 1):
#>           marker_1  marker_2 marker_3 marker_4
#> marker_1 3715.1950 -280.3326   9.0840   2.0554
#> marker_2 -280.3326 1789.9094 -48.7702   1.1472
#> marker_3    9.0840  -48.7702  91.9307   0.3765
#> marker_4    2.0554    1.1472   0.3765   0.0618
#> 
#> Covariance (Class 2):
#>          marker_1 marker_2 marker_3 marker_4
#> marker_1 727.8720  -1.4997  -1.6707   0.0540
#> marker_2  -1.4997 548.9474   1.0490  -0.0899
#> marker_3  -1.6707   1.0490  28.2157  -0.0230
#> marker_4   0.0540  -0.0899  -0.0230   0.0098
#> 
#> Concomitant Coefficients (reference = Class 1):
#> 
#>   -> Class 2:
#>              Estimate        SE      z        p    
#> (Intercept)  3.490840  0.455513  7.664 1.81e-14 ***
#> age         -0.047073  0.009305 -5.059 4.21e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
# }
```
