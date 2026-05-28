# Classification Uncertainty Histogram

Displays the distribution of maximum posterior probabilities across
observations. A peak near 1.0 indicates confident classification; mass
near 1/K indicates substantial ambiguity.

## Usage

``` r
.plot_uncertainty(model, ...)
```

## Arguments

- model:

  A `mixLCA` object.

- ...:

  Additional arguments passed to
  [`geom_histogram`](https://ggplot2.tidyverse.org/reference/geom_histogram.html)
  (e.g. `bins`, `binwidth`).

## Value

A `ggplot` object.
