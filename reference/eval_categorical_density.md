# Evaluate Categorical Log-Density

Computes the log-probability of each observation's categorical
indicators under class-specific multinomial parameters. Missing values
are marginalised out (contribute zero to the sum).

## Usage

``` r
eval_categorical_density(df, probs)
```

## Arguments

- df:

  Data frame of categorical indicators (N x J).

- probs:

  Named list: each element is a named probability vector for one
  variable.

## Value

Numeric vector of length N (log-densities).
