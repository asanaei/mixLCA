# Loadings Plot for a Spectral Dimension

Loadings Plot for a Spectral Dimension

## Usage

``` r
.plot_spectral_loadings(model, dimension = 1L, class = NULL, n_top = NULL)
```

## Arguments

- model:

  A `mixLCA` object.

- dimension:

  Integer in `1:d`.

- class:

  Integer class index.

- n_top:

  Optional integer: keep only top absolute loadings.

## Value

A `ggplot` object.
