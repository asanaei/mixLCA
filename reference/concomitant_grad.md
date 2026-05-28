# Analytical Gradient for Concomitant Model

Analytical Gradient for Concomitant Model

## Usage

``` r
concomitant_grad(par, X, posteriors)
```

## Arguments

- par:

  Numeric vector of length P\*(K-1).

- X:

  Numeric matrix (N x P).

- posteriors:

  Numeric matrix (N x K).

## Value

Numeric vector (gradient of the negative log-likelihood).
