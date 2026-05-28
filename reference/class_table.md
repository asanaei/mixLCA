# Classification Table for a mixLCA Model

Cross-tabulates modal class assignments with average posterior
probabilities. Rows = modal assignment, columns = average P(class k).
Each row sums to 1. Diagonal dominance indicates good class separation.

## Usage

``` r
class_table(model)
```

## Arguments

- model:

  A `mixLCA` object.

## Value

K x K matrix: rows = modal assignment, columns = average P(class k).
