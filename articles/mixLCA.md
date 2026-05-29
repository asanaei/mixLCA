# Introduction to mixLCA

## What `mixLCA` does

`mixLCA` fits finite mixture models for latent class analysis with a
**partitioned architecture**: concomitant predictors enter through a
multinomial logistic regression on class membership; manifest indicators
(continuous, categorical, or mixed) define the measurement model; distal
outcomes are estimated *after* the measurement model converges, under
Bolck-Croon-Hagenaars inverse-classification-error weighting. The
partition prevents distal outcomes from leaking into class meaning and
concomitants from contaminating the measurement structure.

The package supports two responses to local dependence among manifest
indicators:

- **Direct effects** for categorical items (Hagenaars-Vermunt
  specification search via
  [`auto_bvr()`](https://asanaei.github.io/mixLCA/reference/auto_bvr.md)).
- **Spectral Local Dependence (SLD)** for categorical items: a low-rank
  decomposition of the conditional Burt matrix
  ([`auto_sld()`](https://asanaei.github.io/mixLCA/reference/auto_sld.md)).

For continuous indicators, the covariance structure can be fully
unrestricted, diagonal, or sparsified via graphical lasso
(`dependence = "penalized"`).

## A minimal example

``` r

library(mixLCA)
data(health_screening)

fit <- fit_lca(health_screening,
               continuous  = c("marker_1", "marker_2", "marker_3", "marker_4"),
               concomitant = ~ age,
               n_classes   = 2,
               control     = lca_control(n_starts = 10))

summary(fit, data = health_screening)
plot(fit, type = "profiles")

distal(fit, health_screening, outcome ~ age, family = "binomial")
```

## Further reading

The package ships three companion vignettes:

- **[`vignette("workflow")`](https://asanaei.github.io/mixLCA/articles/workflow.md)**
  walks through a categorical LCA from scratch on the
  [`poLCA::election`](https://rdrr.io/pkg/poLCA/man/election.html) data:
  naive fit, BVR diagnostics,
  [`auto_bvr()`](https://asanaei.github.io/mixLCA/reference/auto_bvr.md)
  specification search, and SLD comparison. Read this first if your data
  are categorical.

- **[`vignette("distal-covariates")`](https://asanaei.github.io/mixLCA/articles/distal-covariates.md)**
  covers concomitant predictors (character vector vs. formula), the
  [`predict()`](https://rdrr.io/r/stats/predict.html) contract with
  `na.action`-style padding, and BCH distal estimation on the Pima
  diabetes data.

- **[`vignette("sld-theory")`](https://asanaei.github.io/mixLCA/articles/sld-theory.md)**
  explains the mathematics of Spectral Local Dependence: the conditional
  Burt matrix, the rank-$`d`$ projector, identification, and the M-step.
  Read this if you want to understand what SLD is doing and when it is
  the right tool.

## Function index

Core fitting:

- [`fit_lca()`](https://asanaei.github.io/mixLCA/reference/fit_lca.md),
  [`lca_control()`](https://asanaei.github.io/mixLCA/reference/lca_control.md):
  entry point and optimizer settings.
- [`auto_bvr()`](https://asanaei.github.io/mixLCA/reference/auto_bvr.md):
  BVR-guided specification search.
- [`auto_sld()`](https://asanaei.github.io/mixLCA/reference/auto_sld.md):
  adaptive Spectral Local Dependence rank selection.
- [`enumerate_lca()`](https://asanaei.github.io/mixLCA/reference/enumerate_lca.md):
  fit a range of $`K`$ values.

Diagnostics:

- [`fit_indices()`](https://asanaei.github.io/mixLCA/reference/fit_indices.md):
  AIC, BIC, aBIC, entropy, ICL.
- [`bvr_tests()`](https://asanaei.github.io/mixLCA/reference/bvr_tests.md),
  [`bvr_categorical()`](https://asanaei.github.io/mixLCA/reference/bvr_categorical.md):
  bivariate residual diagnostics.
- [`class_table()`](https://asanaei.github.io/mixLCA/reference/class_table.md):
  modal-class x posterior cross-tabulation.

Inference:

- [`predict()`](https://rdrr.io/r/stats/predict.html) (S3 method):
  posteriors, modal class, full diagnostic frame.
- [`distal()`](https://asanaei.github.io/mixLCA/reference/distal.md):
  BCH-weighted regression for distal outcomes.
- [`concomitant_se()`](https://asanaei.github.io/mixLCA/reference/concomitant_se.md),
  [`continuous_se()`](https://asanaei.github.io/mixLCA/reference/continuous_se.md):
  asymptotic standard errors.
- [`get_posteriors()`](https://asanaei.github.io/mixLCA/reference/get_posteriors.md),
  [`get_loadings()`](https://asanaei.github.io/mixLCA/reference/get_loadings.md),
  [`coef()`](https://rdrr.io/r/stats/coef.html): accessors.

Visualization:

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) (S3 method)
  with
  `type = c("profiles", "bvr", "distal", "uncertainty", "convergence", "categorical", "spectral_scree", "spectral_loadings")`.

## Session info

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.57         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.0     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.6.0    tools_4.6.0       ragg_1.5.2        evaluate_1.0.5   
#> [21] bslib_0.11.0      yaml_2.3.12       jsonlite_2.0.0    rlang_1.2.0      
#> [25] fs_2.1.0
```
