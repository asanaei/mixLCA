# EM Engine for mixLCA

Runs the full EM algorithm: initialises parameters via k-means, iterates
E/M steps until convergence, and returns a fitted `mixLCA` object.

## Usage

``` r
run_em(
  data,
  continuous = NULL,
  categorical = NULL,
  concomitant = NULL,
  n_classes = 2L,
  dependence = "full",
  penalty = 0,
  max_iter = 500L,
  tol = 1e-06,
  seed = 110L,
  cat_direct_effects = NULL,
  spectral_rank = 0L,
  spectral_pool = FALSE,
  init_model = NULL,
  kmeans_nstart = 1L
)
```

## Arguments

- data:

  Data frame.

- continuous:

  Character vector of continuous indicator names, or NULL.

- categorical:

  Character vector of categorical indicator names, or NULL.

- concomitant:

  Character vector of concomitant predictor names, or NULL.

- n_classes:

  Integer \>= 2.

- dependence:

  One of `"none"`, `"full"`, `"penalized"`.

- penalty:

  Numeric L1 penalty for penalised dependence.

- max_iter:

  Maximum EM iterations.

- tol:

  Convergence tolerance on absolute log-likelihood change.

- seed:

  Random seed for initialisation.

- cat_direct_effects:

  List of direct effect pairs, or NULL.

- spectral_rank:

  Integer: target rank for SLD.

- spectral_pool:

  Logical: pool Burt matrices across classes.

- init_model:

  Optional prior `mixLCA` object for warm-start.

## Value

An object of class `mixLCA`.
