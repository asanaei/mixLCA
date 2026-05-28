# Compute BCH Observation Weights

Inverts the classification error matrix and produces per-observation,
per-class weights that correct for classification uncertainty.

## Usage

``` r
bch_weights(posteriors)
```

## Arguments

- posteriors:

  N x K posterior probability matrix.

## Value

List with elements `weights` (N x K), `W_inv` (K x K inverse error
matrix), and `error_matrix` (K x K).
