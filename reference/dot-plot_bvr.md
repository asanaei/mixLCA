# Bivariate Residual Network Graph

Constructs a network graph where nodes are continuous indicators and
edges connect pairs whose residual covariance is significant at p \< .05
(default) or exceeds a user-supplied numeric threshold. Edge thickness
and opacity encode residual magnitude; edge colour distinguishes
positive (red) from negative (blue) residuals.

## Usage

``` r
.plot_bvr(model, data, threshold = "sig")
```

## Arguments

- model:

  A `mixLCA` object.

- data:

  Data frame used for estimation.

- threshold:

  Numeric minimum absolute residual covariance, or `"sig"` to filter by
  chi-squared p \< .05.

## Value

A `ggplot` object via `ggraph`, or NULL when no pairs exceed the
threshold.
