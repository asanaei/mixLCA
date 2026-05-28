# Package index

## Package overview

- [`mixLCA`](https://asanaei.github.io/mixLCA/reference/mixLCA-package.md)
  [`mixLCA-package`](https://asanaei.github.io/mixLCA/reference/mixLCA-package.md)
  : mixLCA: Partitioned Latent Class Analysis

## Fitting

Entry points for fitting latent class models.

- [`fit_lca()`](https://asanaei.github.io/mixLCA/reference/fit_lca.md) :
  Fit a Latent Class Model with mixLCA
- [`lca_control()`](https://asanaei.github.io/mixLCA/reference/lca_control.md)
  : Control Parameters for mixLCA
- [`enumerate_lca()`](https://asanaei.github.io/mixLCA/reference/enumerate_lca.md)
  : Enumerate mixLCA Models Across Class Counts
- [`auto_bvr()`](https://asanaei.github.io/mixLCA/reference/auto_bvr.md)
  : Automated Model Selection for mixLCA
- [`auto_sld()`](https://asanaei.github.io/mixLCA/reference/auto_sld.md)
  : Automated Class-Specific Spectral Rank Selection

## Inference

Post-fit inference and prediction.

- [`predict(`*`<mixLCA>`*`)`](https://asanaei.github.io/mixLCA/reference/predict.mixLCA.md)
  : Predict Method for mixLCA
- [`distal()`](https://asanaei.github.io/mixLCA/reference/distal.md) :
  Distal Outcome Estimation via BCH Weighting
- [`concomitant_se()`](https://asanaei.github.io/mixLCA/reference/concomitant_se.md)
  : Standard Errors for Concomitant Coefficients
- [`continuous_se()`](https://asanaei.github.io/mixLCA/reference/continuous_se.md)
  : Standard Errors for Continuous Parameters (Numerical)
- [`get_posteriors()`](https://asanaei.github.io/mixLCA/reference/get_posteriors.md)
  : Extract Posterior Class Probabilities
- [`get_loadings()`](https://asanaei.github.io/mixLCA/reference/get_loadings.md)
  : Extract Spectral Loadings
- [`coef(`*`<mixLCA>`*`)`](https://asanaei.github.io/mixLCA/reference/coef.mixLCA.md)
  : Coef Method for mixLCA
- [`logLik(`*`<mixLCA>`*`)`](https://asanaei.github.io/mixLCA/reference/logLik.mixLCA.md)
  : logLik Method for mixLCA
- [`nobs(`*`<mixLCA>`*`)`](https://asanaei.github.io/mixLCA/reference/nobs.mixLCA.md)
  : nobs Method for mixLCA
- [`AIC(`*`<mixLCA>`*`)`](https://asanaei.github.io/mixLCA/reference/AIC.mixLCA.md)
  : AIC Method for mixLCA
- [`BIC(`*`<mixLCA>`*`)`](https://asanaei.github.io/mixLCA/reference/BIC.mixLCA.md)
  : BIC Method for mixLCA

## Diagnostics

Fit indices and bivariate residual diagnostics.

- [`fit_indices()`](https://asanaei.github.io/mixLCA/reference/fit_indices.md)
  : Compute Fit Indices for a mixLCA Model
- [`bvr()`](https://asanaei.github.io/mixLCA/reference/bvr.md) :
  Bivariate Residual Covariance Matrix
- [`bvr_tests()`](https://asanaei.github.io/mixLCA/reference/bvr_tests.md)
  : Bivariate Residual Significance Tests
- [`bvr_categorical()`](https://asanaei.github.io/mixLCA/reference/bvr_categorical.md)
  : Bivariate Residuals for Categorical Indicators
- [`class_table()`](https://asanaei.github.io/mixLCA/reference/class_table.md)
  : Classification Table for a mixLCA Model
- [`compare_models()`](https://asanaei.github.io/mixLCA/reference/compare_models.md)
  : Compare Multiple mixLCA Models
- [`spectral_loadings()`](https://asanaei.github.io/mixLCA/reference/spectral_loadings.md)
  : Spectral Loadings as a Tidy Data Frame

## Visualization

Plot method and visual diagnostics.

- [`plot(`*`<mixLCA>`*`)`](https://asanaei.github.io/mixLCA/reference/plot.mixLCA.md)
  : Plot Method for mixLCA

## Printing

- [`print(`*`<mixLCA>`*`)`](https://asanaei.github.io/mixLCA/reference/print.mixLCA.md)
  : Print a mixLCA Object
- [`summary(`*`<mixLCA>`*`)`](https://asanaei.github.io/mixLCA/reference/summary.mixLCA.md)
  : Summary Method for mixLCA
