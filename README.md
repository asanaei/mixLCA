# mixLCA

Partitioned latent class analysis for R, with mixed-type indicators,
concomitant predictors, distal outcomes, and a low-rank treatment of
local dependence (Spectral Local Dependence, SLD).

Website: <https://asanaei.github.io/mixLCA/>

## What it does

- **Mixed indicators.** Continuous (multivariate normal with optional
  graphical-lasso penalisation) and categorical (product-multinomial),
  jointly or alone.
- **Partitioned architecture.** Antecedent concomitant predictors enter
  via multinomial logistic regression on class membership;
  contemporaneous manifest indicators define the measurement model;
  subsequent distal outcomes are estimated separately under
  Bollen-Croon-Hagenaars (BCH) inverse-classification-error weighting.
- **Local dependence.** Direct effects between categorical pairs
  (Hagenaars-Vermunt specification search via `auto_bvr()`) and
  Spectral Local Dependence (`auto_sld()`), a rank-d eigendecomposition
  of the conditional Burt matrix.
- **Honest standard errors.** Full-Hessian SEs for continuous
  parameters, sandwich SEs with NA-on-failure for distal regressions.

## Installation

```r
# install.packages("remotes")
remotes::install_github("asanaei/mixLCA")
```

The package compiles C++ via Rcpp and RcppArmadillo, so you will need a
working compiler toolchain.

## Quick start

```r
library(mixLCA)

data("PimaIndiansDiabetes2", package = "mlbench")
pima <- na.omit(PimaIndiansDiabetes2)

fit <- fit_lca(pima,
               continuous  = c("glucose", "pressure", "mass", "pedigree"),
               concomitant = ~ age,
               n_classes   = 2,
               control     = lca_control(n_starts = 10, seed = 110))

summary(fit, data = pima)
plot(fit, type = "profiles")

distal(fit, pima, diabetes ~ age, family = "binomial")
```

## Vignettes

The website hosts three companion vignettes:

- **A workflow for LCA**: naive fit → BVR diagnostics → direct effects → SLD.
- **Concomitant predictors and distal outcomes**: covariates (character vector and formula), `predict()` semantics, BCH distal estimation.
- **SLD theory and worked example**: the math (conditional Burt matrix, hollow projector, identification, EM step) plus a synthetic example where SLD outperforms naive LCA.

## Citation

The methodological paper for SLD is in preparation. For now, please
cite the package:

```r
citation("mixLCA")
```

## Author

Ali Sanaei (`sanaei@uchicago.edu`).
