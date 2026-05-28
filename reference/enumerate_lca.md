# Enumerate mixLCA Models Across Class Counts

Fits
[`fit_lca()`](https://asanaei.github.io/mixLCA/reference/fit_lca.md) for
each value in `k_range` and returns all fitted objects together with a
model comparison table. To run the starts in parallel, set a
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
in your session before calling this function.

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

  Character vector or formula or NULL.

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

## Examples

``` r
if (FALSE) { # \dontrun{
data(voter_perceptions)
enum <- enumerate_lca(
  data        = voter_perceptions,
  categorical = names(voter_perceptions),
  k_range     = 2:4,
  n_starts    = 5,
  verbose     = FALSE)
enum$comparison
} # }
```
