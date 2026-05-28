# Bivariate Residuals for Categorical Indicators

Computes the Pearson chi-squared bivariate residual for each pair of
categorical manifest variables. Under local independence the
model-implied bivariate frequencies should match the observed ones.
Large BVR values (rule of thumb \> 4) signal local dependence
addressable by adding direct effects via `cat_direct_effects` in `lca`.

## Usage

``` r
bvr_categorical(model, data)
```

## Arguments

- model:

  A `mixLCA` object.

- data:

  Data frame used for estimation.

## Value

Data frame with columns `var1`, `var2`, `bvr`, `df`, `p_value`, ordered
by descending BVR.
