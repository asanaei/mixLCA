# Soft-Threshold Operator

Applies element-wise soft thresholding: sign(x) \* max(\|x\| - lambda,
0). Used for L1-penalised covariance estimation in mixLCA.

## Usage

``` r
soft_threshold(x, lambda)
```

## Arguments

- x:

  Numeric value or matrix.

- lambda:

  Non-negative penalty.

## Value

Thresholded value(s).
