# Multi-Start EM Estimation for mixLCA

Runs the EM algorithm from multiple random starting points (each seeded
deterministically from `base_seed + s`) and returns the solution with
the highest terminal log-likelihood.

## Usage

``` r
run_em_robust(
  data,
  continuous = NULL,
  categorical = NULL,
  concomitant = NULL,
  n_classes = 2L,
  dependence = "full",
  penalty = 0,
  max_iter = 500L,
  tol = 1e-06,
  n_starts = 1L,
  base_seed = 110L,
  verbose = TRUE,
  n_cores = 1L,
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

  Character vector or NULL.

- categorical:

  Character vector or NULL.

- concomitant:

  Character vector or NULL.

- n_classes:

  Integer \>= 2.

- dependence:

  Character.

- penalty:

  Numeric.

- max_iter:

  Integer.

- tol:

  Numeric.

- n_starts:

  Number of random starting configurations.

- base_seed:

  Base random seed.

- verbose:

  Logical: print per-start progress?

- n_cores:

  Integer: number of cores for parallel execution.

- cat_direct_effects:

  List of direct effect pairs, or NULL.

- spectral_rank:

  Integer: SLD rank (0 = disabled).

- spectral_pool:

  Logical: pool Burt matrices across classes.

- init_model:

  Optional prior `mixLCA` object for warm-start (used for the first
  start only).

## Value

The best-fitting `mixLCA` object.
