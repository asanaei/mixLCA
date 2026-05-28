# Control Parameters for mixLCA

Bundles optimiser settings into a single list for use with
[`fit_lca`](https://asanaei.github.io/mixLCA/reference/fit_lca.md).

## Usage

``` r
lca_control(max_iter = 500L, tol = 1e-6, n_starts = 1L,
            seed = 110L, kmeans_nstart = 1L)
```

## Arguments

- max_iter:

  Maximum EM iterations per start.

- tol:

  Convergence tolerance on absolute log-likelihood change.

- n_starts:

  Integer: number of random starting points.

- seed:

  Base random seed; start *s* uses `seed + s`. The global `.Random.seed`
  is never modified.

- kmeans_nstart:

  Integer: random initializations for the internal `kmeans` used to seed
  continuous-indicator starting values.

## Value

A list of control values to pass to
[`fit_lca`](https://asanaei.github.io/mixLCA/reference/fit_lca.md).
