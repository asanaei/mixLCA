# Latent Class Clustering

`lca_clust()` defines a latent class model for clustering via finite
mixture estimation.

This specification is designed for use with the
[tidyclust](https://tidyclust.tidymodels.org) framework. The
engine-specific details for the mixLCA engine are documented below.

## Usage

``` r
lca_clust(
  mode = "partition",
  engine = "mixLCA",
  num_clusters = NULL,
  dependence = NULL,
  penalty = NULL,
  spectral_rank = NULL
)
```

## Arguments

- mode:

  A single character string for the type of model. The only possible
  value for this model is `"partition"`.

- engine:

  A single character string specifying the computational engine.
  Currently only `"mixLCA"`.

- num_clusters:

  Positive integer, the number of latent classes (required).

- dependence:

  Character controlling the covariance structure of continuous
  indicators: `"none"` (diagonal), `"full"` (unrestricted), or
  `"penalized"` (graphical-lasso).

- penalty:

  Non-negative numeric penalty for graphical-lasso estimation (used when
  `dependence = "penalized"`), or `"auto"`. Can be set to `tune()` for
  grid search.

- spectral_rank:

  Non-negative integer rank for Spectral Local Dependence. Zero means
  standard (no SLD). Can be set to `tune()` for grid search.

## Value

An `lca_clust` cluster specification.

## Details

\## What does it mean to predict?

To predict the cluster assignment for a new observation, the model
computes posterior class probabilities under the estimated mixture and
returns the modal (most probable) class.

\## Engine arguments

Additional arguments may be passed to the underlying
[`fit_lca()`](https://asanaei.github.io/mixLCA/reference/fit_lca.md)
function via `set_engine("mixLCA", ...)`. Commonly used engine arguments
include:

- `continuous`:

  Character vector of continuous indicator column names.

- `categorical`:

  Character vector of categorical indicator column names.

- `concomitant`:

  Character vector or one-sided formula for concomitant predictors.

- `penalty`:

  Non-negative numeric or `"auto"`.

- `spectral_rank`:

  Integer for Spectral Local Dependence rank.

- `control`:

  A list from
  [`lca_control()`](https://asanaei.github.io/mixLCA/reference/lca_control.md).

When neither `continuous` nor `categorical` is specified, column types
are inferred automatically: numeric columns become continuous indicators
and factor or character columns become categorical indicators.

## Examples

``` r
# \donttest{
if (requireNamespace("tidyclust", quietly = TRUE)) {
  data(voter_perceptions)
  set.seed(110)

  # --- Naive LCA (local independence) ---
  spec_naive <- lca_clust(num_clusters = 3) |>
    tidyclust::set_engine("mixLCA",
      categorical = names(voter_perceptions),
      control = lca_control(n_starts = 3))

  fit_naive <- tidyclust::fit(spec_naive, ~ ., data = voter_perceptions)
  predict(fit_naive, new_data = voter_perceptions)

  # --- SLD clustering (spectral local dependence, rank 2) ---
  spec_sld <- lca_clust(num_clusters = 3, spectral_rank = 2L) |>
    tidyclust::set_engine("mixLCA",
      categorical = names(voter_perceptions),
      control = lca_control(n_starts = 3))

  fit_sld <- tidyclust::fit(spec_sld, ~ ., data = voter_perceptions)
  predict(fit_sld, new_data = voter_perceptions)

  # Compare fit indices
  cat("Naive BIC:", BIC(fit_naive$fit),
      " SLD BIC:", BIC(fit_sld$fit), "\n")
}
#> Note: `dependence` and `penalty` only affect continuous indicators. With categorical-only models they have no effect.
#> Note: `dependence` and `penalty` only affect continuous indicators. With categorical-only models they have no effect.
#> Naive BIC: 61256.36  SLD BIC: 61791.25 
# }
```
