# ==============================================================================
# File: R/mixLCA-package.R
# Package-level documentation and import declarations for roxygen2.
# ==============================================================================

#' mixLCA: Partitioned Latent Class Analysis
#'
#' Finite mixture modelling with partitioned architecture for antecedent
#' predictors, contemporaneous manifest indicators, and distal outcomes.
#'
#' @useDynLib mixLCA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats AIC BIC as.formula cov kmeans median model.frame
#'   model.matrix model.response na.omit na.pass optim pchisq pnorm
#'   printCoefmat rnorm runif var complete.cases binomial poisson
#'   logLik nobs
#' @importFrom ggplot2 aes geom_col geom_density geom_errorbar
#'   geom_histogram geom_line geom_point geom_vline ggplot labs
#'   theme_minimal theme_void
"_PACKAGE"
