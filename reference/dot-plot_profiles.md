# Profile Plot of Continuous Indicator Means

Renders class-specific mean profiles with optional 95% confidence
intervals derived from numerical standard errors. Classes are optionally
reordered by profile severity (sum of means) to produce a consistent
severity gradient across model runs.

## Usage

``` r
.plot_profiles(model, data = NULL, ci = FALSE, reorder = TRUE, ...)
```

## Arguments

- model:

  A `mixLCA` object with continuous indicators.

- data:

  Data frame (required if `ci = TRUE`).

- ci:

  Logical: draw 95% confidence intervals?

- reorder:

  Logical: reorder classes by profile severity?

- ...:

  Additional arguments passed to
  [`geom_line`](https://ggplot2.tidyverse.org/reference/geom_path.html)
  (e.g. `linewidth`).

## Value

A `ggplot` object.
