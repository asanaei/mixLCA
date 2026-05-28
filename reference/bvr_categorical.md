# Bivariate Residuals for Categorical Indicators

Computes the Pearson chi-squared bivariate residual for each pair of
categorical manifest variables. Under local independence the
model-implied bivariate frequencies should match the observed ones.
Large BVR values (rule of thumb \> 4) signal local dependence
addressable by adding direct effects via `cat_direct_effects` in
[`fit_lca`](https://asanaei.github.io/mixLCA/reference/fit_lca.md).

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

## Examples

``` r
# \donttest{
data(voter_perceptions)
fit <- fit_lca(voter_perceptions,
               categorical = names(voter_perceptions),
               n_classes   = 3,
               control     = lca_control(n_starts = 2, seed = 110),
               verbose     = FALSE)
head(bvr_categorical(fit, voter_perceptions), 5)
#>        var1     var2      bvr df      p_value
#> 60 honest_B empath_B 39.06985  9 1.119123e-05
#> 55 honest_A empath_A 36.92612  9 2.712931e-05
#> 10  moral_A empath_A 36.18033  9 3.682914e-05
#> 53 honest_A  intel_A 34.90915  9 6.182226e-05
#> 6   moral_A honest_A 32.06818  9 1.937961e-04
# }
```
