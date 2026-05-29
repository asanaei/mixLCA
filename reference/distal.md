# Distal Outcome Estimation via BCH Weighting

Estimates class-specific regression models for one or more distal
outcomes, using BCH inverse-classification-error weights. The
measurement model posteriors are fixed before this step (the "cut"): no
gradient from distal outcomes reaches the class definitions.

## Usage

``` r
distal(model, data, formula, family = "gaussian")
```

## Arguments

- model:

  Fitted `mixLCA` object.

- data:

  Data frame containing both the original variables and the distal
  outcome.

- formula:

  A formula for the distal model, e.g.  
  `outcome ~ predictor1 + predictor2`.  
  A right-hand side of `~ 1` estimates unconditional class means.

- family:

  Character: `"gaussian"`, `"binomial"`, or `"poisson"`.

## Value

An object of class `mixDistal` containing class-specific model
summaries.

## Examples

``` r
# \donttest{
data(health_screening)
fit <- fit_lca(health_screening,
               continuous  = c("marker_1","marker_2","marker_3","marker_4"),
               concomitant = ~ age,
               n_classes   = 2,
               control     = lca_control(n_starts = 3, seed = 110),
               verbose     = FALSE)

# Binary distal outcome under BCH weighting
d_bin <- distal(fit, health_screening,
                outcome ~ age, family = "binomial")
print(d_bin)
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

# Continuous distal outcome (uses marker_4 as a stand-in)
d_gauss <- distal(fit, health_screening,
                  marker_4 ~ age, family = "gaussian")
print(d_gauss)
#> 
#> Distal Outcome Estimation (BCH Method) - mixLCA
#> ================================================
#> Formula: marker_4 ~ age 
#> Family : gaussian 
#> Classes: 2 
#> 
#> --- Class 1 ---
#>               Estimate         SE      z        p    
#> (Intercept)  0.6334106  0.1072167  5.908 3.47e-09 ***
#> age         -0.0008039  0.0020789 -0.387    0.699    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>   Residual SD: 0.2432 
#>   Effective N: 167.3 
#> 
#> --- Class 2 ---
#>              Estimate        SE      z      p    
#> (Intercept)  0.359361  0.018619 19.301 <2e-16 ***
#> age         -0.001429  0.000438 -3.264 0.0011 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>   Residual SD: 0.0544 
#>   Effective N: 632.7 
#> 
#> Classification Error Matrix:
#>        [,1]  [,2]
#> [1,] 0.8134 0.022
#> [2,] 0.1866 0.978
#> 
# }
```
