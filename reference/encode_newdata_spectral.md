# Re-Encode New Data Using Training Spectral Encoding

Builds the indicator matrix Z and missingness mask for out-of-sample
data, using the category levels fixed at training time. Unseen
categories are treated as missing.

## Usage

``` r
encode_newdata_spectral(df_new, encoding)
```

## Arguments

- df_new:

  Data frame of categorical indicators.

- encoding:

  Training-time encoding list from `model$cat_spectral_params$encoding`.

## Value

List with elements `Z` and `Z_mis`.
