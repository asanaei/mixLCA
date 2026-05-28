# Validate Inputs for lca

Checks data, variable names, and parameter sanity before estimation.

## Usage

``` r
validate_inputs(
  data,
  continuous,
  categorical,
  concomitant,
  n_classes,
  dependence,
  penalty
)
```

## Arguments

- data:

  Data frame.

- continuous:

  Character vector or NULL.

- categorical:

  Character vector or NULL.

- concomitant:

  Character vector or NULL.

- n_classes:

  Positive integer.

- dependence:

  One of `"none"`, `"full"`, `"penalized"`.

- penalty:

  Non-negative numeric.

## Value

Invisible TRUE on success; stops with an informative message otherwise.
