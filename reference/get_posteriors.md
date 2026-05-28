# Extract Posterior Class Probabilities

Accessor returning the posterior probability matrix. Use this rather
than `model$posteriors` so downstream code is insulated against future
restructuring of internal fields.

## Usage

``` r
get_posteriors(x, ...)

# S3 method for class 'mixLCA'
get_posteriors(x, ...)
```

## Arguments

- x:

  A fitted `mixLCA` object.

- ...:

  Unused.

## Value

Numeric matrix of posterior class probabilities (N x K).
