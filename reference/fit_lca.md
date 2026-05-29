# Fit a Latent Class Model with mixLCA

Estimates a finite mixture model with a partitioned architecture.
Antecedent concomitant predictors enter via multinomial logistic
regression; contemporaneous manifest indicators (continuous,
categorical, or mixed) define the measurement model; subsequent distal
outcomes are estimated separately via
[`distal`](https://asanaei.github.io/mixLCA/reference/distal.md) under
BCH inverse-classification-error weighting.

## Usage

``` r
fit_lca(
  data,
  continuous = NULL,
  categorical = NULL,
  concomitant = NULL,
  n_classes = 2L,
  dependence = "full",
  penalty = "auto",
  control = lca_control(),
  verbose = TRUE,
  cat_direct_effects = NULL,
  spectral_rank = 0L,
  spectral_pool = FALSE,
  init_model = NULL
)
```

## Arguments

- data:

  Data frame of observations.

- continuous:

  Character vector of continuous indicator variable names, or NULL.

- categorical:

  Character vector of categorical indicator variable names, or NULL.

- concomitant:

  Character vector of concomitant predictor names, or a one-sided
  formula (e.g.\\ `~ age * income + poly(bmi, 2)`), or NULL. Missing
  values in concomitant predictors are not permitted; impute or filter
  rows before calling.

- n_classes:

  Integer \>= 2: number of latent classes.

- dependence:

  Character controlling the covariance structure of continuous
  indicators within each class:

  `"none"`

  :   Diagonal covariance (local independence).

  `"full"`

  :   Unrestricted covariance.

  `"penalized"`

  :   Graphical-lasso penalized covariance, guaranteeing exact sparsity
      and positive definiteness.

- penalty:

  Penalty for `dependence = "penalized"`. Either a non-negative scalar
  or the string `"auto"` (default). When `"auto"`, a heuristic value is
  selected from the data; when numeric, the supplied value is respected
  exactly (including 0, which then yields no shrinkage).

- control:

  List of optimizer settings; see
  [`lca_control`](https://asanaei.github.io/mixLCA/reference/lca_control.md).

- verbose:

  Logical: print per-start progress?

- cat_direct_effects:

  List of two-element character vectors specifying direct effects
  between categorical indicators to address local dependence violations.
  Each element `c("parent", "child")` allows the child's response
  probabilities to depend on the parent's observed value within each
  class (Vermunt, 1999). Use
  [`bvr_categorical`](https://asanaei.github.io/mixLCA/reference/bvr_categorical.md)
  to identify candidate pairs. When `NULL` (default), standard local
  independence is assumed for categorical indicators.

- spectral_rank:

  Integer (scalar or length-K vector): target rank *d* for Spectral
  Local Dependence (SLD) among categorical items. If any element is \>
  0, SLD is activated.

- spectral_pool:

  Logical: Pool the conditional Burt matrices across classes.

- init_model:

  Optional `mixLCA` object for warm-starting. The prior model must have
  the same `n_classes`.

## Value

An object of class `mixLCA` (a list) containing:

- n_classes:

  The number of latent classes.

- log_lik:

  The maximized log-likelihood at the final EM iterate.

- posteriors:

  An N x K matrix of posterior class probabilities.

- continuous_params:

  List of class-specific means and covariances (when continuous
  indicators are used).

- categorical_params:

  List of class-specific category probabilities (when categorical
  indicators are used).

- cat_spectral_params:

  List of class-specific Spectral Local Dependence parameters (when
  `spectral_rank > 0`).

- concomitant_coefs:

  Matrix of multinomial logistic regression coefficients (when
  concomitant predictors are supplied).

- n_obs:

  Number of observations used in estimation.

- convergence:

  List detailing EM convergence status and history.

- specs:

  List of model specifications (indicator types, dependence structure,
  etc.).

- AIC, BIC, aBIC, ICL, entropy:

  Cached fit indices.

Treat the return as opaque: use the provided accessors and S3 methods
([`coef()`](https://rdrr.io/r/stats/coef.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`get_posteriors()`](https://asanaei.github.io/mixLCA/reference/get_posteriors.md),
[`get_loadings()`](https://asanaei.github.io/mixLCA/reference/get_loadings.md),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html)) rather than
`$`-indexing internal fields, which may be restructured in future
versions. Key downstream functions:
[`distal`](https://asanaei.github.io/mixLCA/reference/distal.md),
[`fit_indices`](https://asanaei.github.io/mixLCA/reference/fit_indices.md),
[`bvr_tests`](https://asanaei.github.io/mixLCA/reference/bvr_tests.md),
[`enumerate_lca`](https://asanaei.github.io/mixLCA/reference/enumerate_lca.md).

## Details

The measurement likelihood is maximized by EM. Continuous indicators
follow class-specific multivariate normal distributions; missing entries
are marginalized analytically. Categorical indicators follow a
product-multinomial distribution; missing categories are omitted from
the likelihood contribution. Non-diagonal covariance matrices may be
freely estimated (`dependence = "full"`) or estimated with sparsity
through the graphical lasso (`dependence = "penalized"`).

Parallel execution of multiple starts is delegated to the user's active
[`future::plan()`](https://future.futureverse.org/reference/plan.html).
Set the plan in your session
(`future::plan(future::multisession, workers = ...)`) before calling
`fit_lca`; otherwise starts run sequentially.

## Examples

``` r
# \donttest{
# Categorical example: fit a two-class model to voter perceptions.
data(voter_perceptions)
fit_cat <- fit_lca(voter_perceptions,
                   categorical = names(voter_perceptions),
                   n_classes   = 2,
                   control     = lca_control(n_starts = 3),
                   verbose     = FALSE)
fit_cat
#> 
#> Latent Class Model - mixLCA
#> ===========================
#> Classes        : 2 
#> Log-likelihood       : -30684.79 
#> Parameters     : 73 
#> Observations   : 2000 
#> Converged      : TRUE (29 iterations) 
#> 
#> Specification:
#>   Categorical indicators: 12 
#>   Dependence structure  : full 
#> 
#> Class proportions:
#> Class 1 Class 2 
#>  0.5152  0.4848 
#> 

# Mixed example: continuous markers with an age covariate.
data(health_screening)
fit_mix <- fit_lca(health_screening,
                   continuous  = c("marker_1", "marker_2", "marker_3", "marker_4"),
                   concomitant = ~ age,
                   n_classes   = 2,
                   control     = lca_control(n_starts = 5),
                   verbose     = FALSE)
summary(fit_mix, data = health_screening)
#> 
#> Summary - mixLCA
#> ================
#> 
#> Fit Indices:
#>   Log-likelihood     : -9938.795 
#>   AIC               : 19937.59 
#>   BIC               : 20078.13 
#>   Sample-adj. BIC   : 19982.86 
#>   Entropy           : 0.783 
#>   ICL               : 20318.83 
#> 
#> Class Proportions:
#> Class 1 Class 2 
#>  0.2092  0.7908 
#> 
#> Continuous Indicator Means:
#>         marker_1 marker_2 marker_3 marker_4
#> Class 1 130.9492  97.6235  31.1184   0.5409
#> Class 2  82.6076  73.9969  22.3791   0.3120
#> 
#> Covariance (Class 1):
#>           marker_1  marker_2 marker_3 marker_4
#> marker_1 3715.1944 -280.3323   9.0840   2.0554
#> marker_2 -280.3323 1789.9091 -48.7702   1.1472
#> marker_3    9.0840  -48.7702  91.9307   0.3765
#> marker_4    2.0554    1.1472   0.3765   0.0618
#> 
#> Covariance (Class 2):
#>          marker_1 marker_2 marker_3 marker_4
#> marker_1 727.8719  -1.4997  -1.6707   0.0540
#> marker_2  -1.4997 548.9473   1.0490  -0.0899
#> marker_3  -1.6707   1.0490  28.2157  -0.0230
#> marker_4   0.0540  -0.0899  -0.0230   0.0098
#> 
#> Concomitant Coefficients (reference = Class 1):
#> 
#>   -> Class 2:
#>              Estimate        SE      z        p    
#> (Intercept)  3.490835  0.455513  7.664 1.81e-14 ***
#> age         -0.047073  0.009305 -5.059 4.22e-07 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
# }
```
