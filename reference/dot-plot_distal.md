# Distal Outcome Density by Class

Plots kernel density estimates of a distal variable split by modal class
assignment. Observations may be weighted by their maximum posterior
probability to reflect classification uncertainty.

## Usage

``` r
.plot_distal(model, data, variable, weighted = TRUE, ...)
```

## Arguments

- model:

  A `mixLCA` object.

- data:

  Data frame containing the distal variable.

- variable:

  Character: name of the distal variable.

- weighted:

  Logical: weight densities by max posterior probability?

- ...:

  Additional arguments passed to
  [`geom_density`](https://ggplot2.tidyverse.org/reference/geom_density.html)
  (e.g. `bw`, `adjust`, `kernel`).

## Value

A `ggplot` object.
