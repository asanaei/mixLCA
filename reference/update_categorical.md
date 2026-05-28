# M-Step Update for Categorical Parameters

Re-estimates class-conditional category probabilities using posterior
weights. Missing values are excluded from both numerator and
denominator.

## Usage

``` r
update_categorical(df, weights)
```

## Arguments

- df:

  Data frame of categorical indicators (N x J).

- weights:

  Numeric vector of length N (posterior weights for this class).

## Value

Named list of named probability vectors (one per variable).
