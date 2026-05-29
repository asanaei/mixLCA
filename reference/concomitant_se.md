# Standard Errors for Concomitant Coefficients

Constructs the empirical (observed) information matrix using the outer
product of per-observation score vectors, then inverts to obtain
asymptotic variances for the multinomial logistic coefficients.

## Usage

``` r
concomitant_se(model, data)
```

## Arguments

- model:

  A `mixLCA` object.

- data:

  Data frame used for estimation.

## Value

A P x (K-1) matrix of standard errors, or NULL if no concomitant
predictors were specified.

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
se <- concomitant_se(fit, health_screening)
round(cbind(Estimate = fit$concomitant_coefs[, 1], SE = se[, 1]), 4)
#>             Estimate     SE
#> (Intercept)  -3.4908 0.4555
#> age           0.0471 0.0093
# }
```
