# Compute BCH Classification Error Matrix

Constructs the K x K matrix E where E\[w, k\] = P(W = w \| C = k), i.e.
the probability that an observation truly belonging to class k is
modally assigned to class w.

## Usage

``` r
classification_error_matrix(posteriors)
```

## Arguments

- posteriors:

  N x K posterior probability matrix.

## Value

K x K classification error matrix.
