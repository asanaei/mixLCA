# Bivariate Residual Covariance Matrix

Computes the residual covariance between continuous indicators after
subtracting the model-implied mixture covariance. Large residuals signal
local dependence violations not captured by the current model. Use
`plot_bvr()` to visualise the residual network.

## Usage

``` r
bvr(model, data)
```

## Arguments

- model:

  A `mixLCA` object.

- data:

  Data frame used for estimation.

## Value

Named matrix of residual covariances among continuous indicators.
