# Enumerate mixLCA Models Across Class Counts

Fits `lca()` for each value in `k_range` and returns all fitted objects
together with a model comparison table.

## Usage

``` r
enumerate_lca(
  data,
  continuous = NULL,
  categorical = NULL,
  concomitant = NULL,
  k_range = 2:4,
  dependence = "full",
  penalty = 0,
  n_starts = 1L,
  max_iter = 500L,
  tol = 1e-06,
  spectral_rank = 0L,
  spectral_pool = FALSE,
  n_cores = 1L,
  verbose = TRUE,
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

- k_range:

  Integer vector of class counts to estimate, e.g.  
  `2:5`.

- dependence:

  Character.

- penalty:

  Numeric.

- n_starts:

  Integer.

- max_iter:

  Integer.

- tol:

  Numeric.

- spectral_rank:

  Integer: SLD rank (0 = disabled).

- spectral_pool:

  Logical: pool Burt matrices across classes.

- n_cores:

  Integer.

- verbose:

  Logical.

- kmeans_nstart:

  Integer: random starts for internal k-means (default 1).

## Value

List with:

- `$models`:

  Named list of `mixLCA` objects.

- `$comparison`:

  Data frame from
  [`compare_models()`](https://asanaei.github.io/mixLCA/reference/compare_models.md).
