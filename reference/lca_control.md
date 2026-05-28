# Control Parameters for mixLCA

Bundles optimiser settings into a single list for use with
[`fit_lca`](https://asanaei.github.io/mixLCA/reference/fit_lca.md).
Defaults reproduce the behaviour of prior beta releases.

## Usage

``` r
lca_control(
  max_iter = 500L,
  tol = 1e-06,
  n_starts = 1L,
  seed = 110L,
  kmeans_nstart = 1L
)
```

## Arguments

- max_iter:

  Maximum EM iterations per start.

- tol:

  Convergence tolerance on absolute log-likelihood change.

- n_starts:

  Integer: number of random starting points. For publication-quality
  results consider at least 10.

- seed:

  Base random seed; start *s* uses `seed + s`. The global `.Random.seed`
  is never modified.

- kmeans_nstart:

  Integer: random initializations for the internal `kmeans` used to seed
  continuous-indicator starting values. Irrelevant for categorical-only
  models or when `init_model` is supplied.

## Value

A list of control values.

## Examples

``` r
# Default control parameters
lca_control()
#> $max_iter
#> [1] 500
#> 
#> $tol
#> [1] 1e-06
#> 
#> $n_starts
#> [1] 1
#> 
#> $seed
#> [1] 110
#> 
#> $kmeans_nstart
#> [1] 1
#> 

# Tighter tolerance, more random starts, deterministic seed
lca_control(max_iter = 1000L, tol = 1e-8,
            n_starts = 10L, seed = 110L)
#> $max_iter
#> [1] 1000
#> 
#> $tol
#> [1] 1e-08
#> 
#> $n_starts
#> [1] 10
#> 
#> $seed
#> [1] 110
#> 
#> $kmeans_nstart
#> [1] 1
#> 
```
