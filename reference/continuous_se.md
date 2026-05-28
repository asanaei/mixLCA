# Standard Errors for Continuous Parameters (Numerical)

Computes approximate standard errors for class-specific means and
diagonal covariance elements via numerical second derivative of the
observed marginal log-likelihood. This respects Louis's Principle by
differentiating through the mixture, not the Q-function.

## Usage

``` r
continuous_se(model, data, step = 1e-04)
```

## Arguments

- model:

  A `mixLCA` object.

- data:

  Data frame used for estimation.

- step:

  Finite difference step size.

## Value

A list with elements `mean_se` (list of K named vectors) and `cov_se`
(list of K d x d matrices), or NULL if no continuous indicators were
specified.
