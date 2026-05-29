# Control Parameters for mixLCA

Bundles optimizer settings into a single list for use with
[`fit_lca`](https://asanaei.github.io/mixLCA/reference/fit_lca.md).
Defaults reproduce the behavior of prior beta releases.

## Usage

``` r
lca_control(max_iter = 500L, tol = 1e-06, n_starts = 1L, kmeans_nstart = 1L)
```

## Arguments

- max_iter:

  Maximum EM iterations per start.

- tol:

  Convergence tolerance on absolute log-likelihood change.

- n_starts:

  Integer: number of random starting points. For publication-quality
  results consider at least 10.

- kmeans_nstart:

  Integer: random initializations for the internal `kmeans` used to seed
  continuous-indicator starting values. Irrelevant for categorical-only
  models or when `init_model` is supplied.

## Value

A list of control values.

## Reproducibility

[`fit_lca()`](https://asanaei.github.io/mixLCA/reference/fit_lca.md)
draws from the global RNG state for random initialization, exactly like
[`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) and
[`uwot::umap()`](https://jlmelville.github.io/uwot/reference/umap.html).
To get a reproducible fit, call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before
[`fit_lca()`](https://asanaei.github.io/mixLCA/reference/fit_lca.md).

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
#> $kmeans_nstart
#> [1] 1
#> 

# Tighter tolerance, more random starts
lca_control(max_iter = 1000L, tol = 1e-8, n_starts = 10L)
#> $max_iter
#> [1] 1000
#> 
#> $tol
#> [1] 1e-08
#> 
#> $n_starts
#> [1] 10
#> 
#> $kmeans_nstart
#> [1] 1
#> 

# Reproducible fit: set seed at the call site, like with kmeans
set.seed(110)
# then call fit_lca(..., control = lca_control(n_starts = 10L))
```
