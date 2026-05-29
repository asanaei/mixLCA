# M-Step: Update Concomitant Coefficients

Wraps BFGS optimization of the multinomial logistic model.

## Usage

``` r
update_concomitant(X, posteriors, current_coefs)
```

## Arguments

- X:

  Numeric matrix (N x P), includes intercept.

- posteriors:

  Numeric matrix (N x K).

- current_coefs:

  Current P x (K-1) coefficient matrix.

## Value

Updated P x (K-1) coefficient matrix.
