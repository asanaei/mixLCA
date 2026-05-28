# Negative Log-Likelihood for Concomitant Model

Negative Log-Likelihood for Concomitant Model

## Usage

``` r
concomitant_nll(par, X, posteriors)
```

## Arguments

- par:

  Numeric vector of length P\*(K-1), stored column-major.

- X:

  Numeric matrix (N x P), includes intercept column.

- posteriors:

  Numeric matrix (N x K) of current posterior probabilities.

## Value

Scalar (negative expected log-likelihood contribution).
