# Automated Model Selection for mixLCA

Three-phase automated model selection:

- Phase 1:

  Select *K* by BIC across `K_range` using diagonal covariance.

- Phase 2a:

  Compare diagonal, penalised, and full covariance structures for
  continuous indicators at the selected *K*.

- Phase 2b:

  Vermunt-style BVR-guided specification search for categorical
  indicators: iteratively add direct effects for the pair with the
  largest bivariate residual, stopping when BIC no longer improves or
  `max_direct_effects` is reached.

## Usage

``` r
auto_bvr(
  data,
  continuous = NULL,
  categorical = NULL,
  concomitant = NULL,
  K_range = 2:5,
  max_direct_effects = 5L,
  bvr_threshold = 3.84,
  seed = 110L,
  n_cores = 4L,
  verbose = TRUE,
  ...
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

- K_range:

  Integer vector of class counts to evaluate (e.g., 2:5).

- max_direct_effects:

  Integer: maximum number of categorical direct effects to add during
  Phase 2b (default 5).

- bvr_threshold:

  Numeric: minimum BVR chi-squared statistic to consider a direct effect
  (default 3.84, i.e. p \< .05 for df = 1).

- seed:

  Random seed.

- n_cores:

  Integer. Number of parallel processes to use. Default is 4.

- verbose:

  Logical: print progress?

- ...:

  Additional arguments passed to `lca`.

## Value

An object of class `mixLCA` representing the final selected model, with
an additional `auto_path` element recording the search trajectory.
