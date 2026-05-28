# EM Convergence Trace Plot

Plots the log-likelihood across EM iterations.

## Usage

``` r
.plot_convergence(model, ...)
```

## Arguments

- model:

  A `mixLCA` object.

- ...:

  Additional arguments passed to
  [`geom_line`](https://ggplot2.tidyverse.org/reference/geom_path.html).

## Value

A `ggplot` object.
