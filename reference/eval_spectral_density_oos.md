# Evaluate SLD Composite Log-Density on New Data

Applies the trained spectral shift to out-of-sample observations.

## Usage

``` r
eval_spectral_density_oos(newdata_encoded, encoding, pi_c, A_star)
```

## Arguments

- newdata_encoded:

  List with Z and Z_mis from `encode_newdata_spectral`.

- encoding:

  Training encoding list.

- pi_c:

  Flattened class-conditional marginal probabilities.

- A_star:

  Hollow projection matrix.

## Value

Numeric vector of length N (log-densities).
