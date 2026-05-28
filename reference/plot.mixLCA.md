# Plot Method for mixLCA

Single dispatch surface for all built-in visualisations.

## Usage

``` r
# S3 method for class 'mixLCA'
plot(x,
                     type = c("profiles", "bvr", "distal",
                              "uncertainty", "convergence",
                              "categorical", "spectral_scree",
                              "spectral_loadings"),
                     ...)
```

## Arguments

- x:

  A fitted `mixLCA` object.

- type:

  One of `"profiles"`, `"bvr"`, `"distal"`, `"uncertainty"`,
  `"convergence"`, `"categorical"`, `"spectral_scree"`, or
  `"spectral_loadings"`.

- ...:

  Forwarded to the underlying plotting routine (e.g.\\ `data` for
  `"bvr"` and `"distal"`, `variable` for `"distal"`, `ci` for
  `"profiles"`, `dimension` for `"spectral_loadings"`).

## Value

A `ggplot` object.
