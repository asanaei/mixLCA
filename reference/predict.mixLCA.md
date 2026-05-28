# Predict Method for mixLCA

Returns posterior class probabilities, modal class assignments, or the
full diagnostic dataframe for each observation. When `newdata` is
supplied, posteriors are computed from scratch using the estimated model
parameters, including Spectral Local Dependence shifts if applicable.

## Usage

``` r
# S3 method for class 'mixLCA'
predict(object, newdata = NULL, type = c("prob", "class", "all"), ...)
```

## Arguments

- object:

  A `mixLCA` object.

- newdata:

  Optional data frame for out-of-sample prediction.

- type:

  One of `"prob"` (default; N x K matrix of posterior probabilities),
  `"class"` (integer vector of modal class assignments), or `"all"`
  (legacy data frame containing `P_class_*`, `modal_class`,
  `max_posterior`, and – for out-of-sample – `log_lik`).

- ...:

  Unused.

## Value

A matrix, integer vector, or data frame, per `type`.

## Details

Rows of `newdata` with missing concomitant values are not scored and are
returned as `NA` so that the output length always equals `nrow(newdata)`
(the standard `na.exclude` contract). This makes
`cbind(newdata, predict(model, newdata))` safe.

## Examples

``` r
# \donttest{
data(voter_perceptions)
fit <- fit_lca(voter_perceptions,
               categorical = names(voter_perceptions),
               n_classes   = 3,
               control     = lca_control(n_starts = 2, seed = 110),
               verbose     = FALSE)

# Default: N x K posterior probability matrix
head(predict(fit))
#>         P_class_1    P_class_2    P_class_3
#> [1,] 1.438034e-05 9.826358e-01 1.734984e-02
#> [2,] 7.110197e-03 7.563908e-06 9.928822e-01
#> [3,] 1.098566e-04 9.992144e-01 6.757417e-04
#> [4,] 1.163068e-05 9.999751e-01 1.322256e-05
#> [5,] 1.748301e-05 9.999775e-01 5.024534e-06
#> [6,] 9.993873e-01 5.790421e-04 3.363004e-05

# Modal class assignments as an integer vector
head(predict(fit, type = "class"))
#> [1] 2 3 2 2 2 1

# Full diagnostic data frame
head(predict(fit, type = "all"))
#>      P_class_1    P_class_2    P_class_3 modal_class max_posterior
#> 1 1.438034e-05 9.826358e-01 1.734984e-02           2     0.9826358
#> 2 7.110197e-03 7.563908e-06 9.928822e-01           3     0.9928822
#> 3 1.098566e-04 9.992144e-01 6.757417e-04           2     0.9992144
#> 4 1.163068e-05 9.999751e-01 1.322256e-05           2     0.9999751
#> 5 1.748301e-05 9.999775e-01 5.024534e-06           2     0.9999775
#> 6 9.993873e-01 5.790421e-04 3.363004e-05           1     0.9993873

# Out-of-sample prediction; rows with NA in concomitants pad with NA
new_rows <- voter_perceptions[1:6, ]
predict(fit, newdata = new_rows, type = "class")
#> [1] 2 3 2 2 2 1
# }
```
