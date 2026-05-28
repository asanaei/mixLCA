# Fit a Latent Class Model with mixLCA

Estimates a finite mixture model with concomitant predictors, mixed-type
measurement indicators, and optional Spectral Local Dependence (SLD).
Parallel multi-start execution honours the user's active
[`future::plan()`](https://future.futureverse.org/reference/plan.html).

## Usage

``` r
fit_lca(
  data,
  continuous = NULL,
  categorical = NULL,
  concomitant = NULL,
  n_classes = 2L,
  dependence = "full",
  penalty = "auto",
  control = lca_control(),
  verbose = TRUE,
  cat_direct_effects = NULL,
  spectral_rank = 0L,
  spectral_pool = FALSE,
  init_model = NULL
)
```

## Arguments

- data:

  Data frame of observations.

- continuous:

  Character vector of continuous indicator variable names, or NULL.

- categorical:

  Character vector of categorical indicator variable names, or NULL.

- concomitant:

  Character vector of concomitant predictor names, or a one-sided
  formula (e.g.\\ `~ age * income + poly(bmi, 2)`), or NULL. Missing
  values in concomitant predictors are not permitted.

- n_classes:

  Integer \>= 2.

- dependence:

  One of `"none"`, `"full"`, or `"penalized"`.

- penalty:

  Non-negative scalar or the string `"auto"`.

- control:

  List from
  [`lca_control`](https://asanaei.github.io/mixLCA/reference/lca_control.md).

- verbose:

  Logical.

- cat_direct_effects:

  List of two-element character vectors specifying parent-\>child direct
  effects between categorical indicators.

- spectral_rank:

  Scalar or length-K integer vector of SLD ranks.

- spectral_pool:

  Logical: pool Burt matrices across classes.

- init_model:

  Optional `mixLCA` object for warm-starting.

## Value

An opaque object of class `mixLCA`.

## See also

[`lca_control`](https://asanaei.github.io/mixLCA/reference/lca_control.md),
[`distal`](https://asanaei.github.io/mixLCA/reference/distal.md),
[`fit_indices`](https://asanaei.github.io/mixLCA/reference/fit_indices.md),
[`bvr_tests`](https://asanaei.github.io/mixLCA/reference/bvr_tests.md),
[`enumerate_lca`](https://asanaei.github.io/mixLCA/reference/enumerate_lca.md).
