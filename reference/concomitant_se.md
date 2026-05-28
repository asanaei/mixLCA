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
