# Plot Categorical Indicator Probabilities

Visualizes the class-conditional response probabilities for categorical
manifest variables using dodged bar charts faceted by variable.

## Usage

``` r
.plot_categorical(
  model,
  variables = NULL,
  orientation = c("vertical", "horizontal")
)
```

## Arguments

- model:

  A `mixLCA` object.

- variables:

  Character vector: optional subset of categorical variables to plot. If
  `NULL`, plots all categorical indicators.

- orientation:

  Character: `"vertical"` (default) or `"horizontal"`. Horizontal is
  recommended for variables with long category names.

## Value

A `ggplot` object.
