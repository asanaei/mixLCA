# Plot Method for mixLCA

Single dispatch surface for all built-in visualisations. The `type`
argument selects which plot to produce.

## Usage

``` r
# S3 method for class 'mixLCA'
plot(
  x,
  y = NULL,
  type = c("profiles", "bvr", "distal", "uncertainty", "convergence", "categorical",
    "spectral_scree", "spectral_loadings"),
  ...
)
```

## Arguments

- x:

  A fitted `mixLCA` object.

- y:

  Unused. Present so the method signature matches the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) generic in
  base R.

- type:

  One of `"profiles"`, `"bvr"`, `"distal"`, `"uncertainty"`,
  `"convergence"`, `"categorical"`, `"spectral_scree"`, or
  `"spectral_loadings"`.

- ...:

  Forwarded to the underlying plotting routine.

## Value

A `ggplot` object.

## Details

Plot-specific arguments (e.g.\\ `data` for `"bvr"` and `"distal"`,
`variable` for `"distal"`, `ci` for `"profiles"`, `dimension` for
`"spectral_loadings"`) are forwarded via `...`.

## Examples

``` r
# \donttest{
# Continuous indicators: profile, uncertainty, convergence
data(health_screening)
fit_c <- fit_lca(health_screening,
                 continuous = c("marker_1","marker_2","marker_3","marker_4"),
                 n_classes  = 2,
                 control    = lca_control(n_starts = 2),
                 verbose    = FALSE)
plot(fit_c, type = "profiles")

plot(fit_c, type = "uncertainty")
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.

plot(fit_c, type = "convergence")

plot(fit_c, type = "bvr", data = health_screening)
#> No bivariate residuals exceed the threshold. Local independence holds adequately.


# Categorical indicators: response probabilities and SLD loadings
data(voter_perceptions)
fit_k <- fit_lca(voter_perceptions,
                 categorical   = names(voter_perceptions)[1:6],
                 n_classes     = 2,
                 spectral_rank = c(1L, 1L),
                 control       = lca_control(n_starts = 2),
                 verbose       = FALSE)
plot(fit_k, type = "categorical")

plot(fit_k, type = "spectral_scree")

plot(fit_k, type = "spectral_loadings", dimension = 1, class = 1)


# Distal density by class
plot(fit_c, type = "distal", data = health_screening,
     variable = "marker_1")

# }
```
