# Extract Spectral Loadings

Accessor returning the SLD loadings table. Use this rather than
`model$cat_spectral_params` so downstream code is insulated against
future restructuring of internal fields.

## Usage

``` r
get_loadings(x, ...)

# S3 method for class 'mixLCA'
get_loadings(x, ...)
```

## Arguments

- x:

  A fitted `mixLCA` object.

- ...:

  Unused.

## Value

A data frame from
[`spectral_loadings`](https://asanaei.github.io/mixLCA/reference/spectral_loadings.md),
or `NULL` if the model has no Spectral Local Dependence component.
