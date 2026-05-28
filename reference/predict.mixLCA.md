# Predict Method for mixLCA

Returns posterior class probabilities, modal class assignments, and the
maximum posterior probability (classification certainty) for each
observation. When `newdata` is supplied, posteriors are computed from
scratch using the estimated model parameters, including Spectral Local
Dependence shifts if applicable.

## Usage

``` r
# S3 method for class 'mixLCA'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  A `mixLCA` object.

- newdata:

  Optional data frame for out-of-sample prediction.

- ...:

  Unused.

## Value

Data frame with one row per observation, including columns `P_class_1`,
..., `P_class_K`, `modal_class`, `max_posterior`, and (when `newdata` is
supplied) `log_lik` (per-observation log marginal density).
