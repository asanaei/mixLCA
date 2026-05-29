# ==============================================================================
# File: R/07_standard_errors.R
# mixLCA - asymptotic standard error estimation.
# ==============================================================================

#' Standard Errors for Concomitant Coefficients
#'
#' Constructs the empirical (observed) information matrix using the
#' outer product of per-observation score vectors, then inverts to
#' obtain asymptotic variances for the multinomial logistic coefficients.
#'
#' @param model A \code{mixLCA} object.
#' @param data Data frame used for estimation.
#' @return A P x (K-1) matrix of standard errors, or NULL if no
#'   concomitant predictors were specified.
#'
#' @examples
#' \donttest{
#' data(health_screening)
#' fit <- fit_lca(health_screening,
#'                continuous  = c("marker_1","marker_2","marker_3","marker_4"),
#'                concomitant = ~ age,
#'                n_classes   = 2,
#'                control     = lca_control(n_starts = 3),
#'                verbose     = FALSE)
#' se <- concomitant_se(fit, health_screening)
#' round(cbind(Estimate = fit$concomitant_coefs[, 1], SE = se[, 1]), 4)
#' }
#' @export
concomitant_se <- function(model, data) {
  if (!inherits(model, "mixLCA"))
    stop("`model` must be a mixLCA object.")
  if (is.null(model$specs$concomitant)) return(NULL)

  coefs      <- model$concomitant_coefs
  posteriors <- model$posteriors
  K  <- model$n_classes
  N  <- model$n_obs

  concom <- model$specs$concomitant
  concom_vars <- if (inherits(concom, "formula")) all.vars(concom) else concom
  concom_complete <- stats::complete.cases(data[, concom_vars, drop = FALSE])
  if (any(!concom_complete)) {
    data <- data[concom_complete, , drop = FALSE]
    posteriors <- posteriors[concom_complete, , drop = FALSE]
    N <- nrow(data)
  }

  if (inherits(concom, "formula")) {
    X <- stats::model.matrix(concom, data = data)
  } else {
    f_str <- paste("~", paste(concom, collapse = " + "))
    X <- stats::model.matrix(stats::as.formula(f_str), data = data)
  }
  P <- ncol(X)

  priors <- compute_priors(X, coefs)

  n_par <- P * (K - 1L)
  S_mat <- matrix(0, nrow = N, ncol = n_par)

  for (i in seq_len(N)) {
    g_i <- matrix(0, nrow = P, ncol = K - 1L)
    for (k in seq_len(K - 1L)) {
      g_i[, k] <- X[i, ] *
        (posteriors[i, k + 1L] - priors[i, k + 1L])
    }
    S_mat[i, ] <- as.vector(g_i)
  }

  I_obs <- crossprod(S_mat)
  eig   <- eigen(I_obs, symmetric = TRUE)
  if (min(eig$values) <= 1e-8)
    I_obs <- I_obs + diag(1e-6, n_par)

  I_inv  <- solve(I_obs)
  se_vec <- sqrt(pmax(diag(I_inv), 0))
  matrix(se_vec, nrow = P, ncol = K - 1L)
}

#' Standard Errors for Continuous Parameters (Numerical)
#'
#' Computes approximate standard errors for class-specific means and
#' diagonal covariance elements via numerical second derivative of the
#' observed marginal log-likelihood. This respects Louis's Principle
#' by differentiating through the mixture, not the Q-function.
#'
#' @param model A \code{mixLCA} object.
#' @param data Data frame used for estimation.
#' @param step Finite difference step size.
#' @return A list with elements \code{mean_se} (list of K named
#'   vectors) and \code{cov_se} (list of K d x d matrices), or NULL
#'   if no continuous indicators were specified.
#'
#' @examples
#' \donttest{
#' data(health_screening)
#' fit <- fit_lca(health_screening,
#'                continuous = c("marker_1","marker_2","marker_3","marker_4"),
#'                n_classes  = 2,
#'                control    = lca_control(n_starts = 3),
#'                verbose    = FALSE)
#' ses <- continuous_se(fit, health_screening)
#' round(ses$mean_se[[1]], 4)   # SE of class-1 means
#' }
#' @export
continuous_se <- function(model, data, step = 1e-4) {
  if (!inherits(model, "mixLCA"))
    stop("`model` must be a mixLCA object.")

  cont <- model$specs$continuous
  if (is.null(cont)) return(NULL)

  Y  <- as.matrix(data[, cont, drop = FALSE])
  K  <- model$n_classes
  d  <- length(cont)

  # Pre-compute categorical density and priors (they stay fixed)
  cat_vars <- model$specs$categorical
  has_cat  <- !is.null(cat_vars)
  D <- if (has_cat) data[, cat_vars, drop = FALSE] else NULL

  has_concom <- !is.null(model$specs$concomitant)
  if (has_concom) {
    concom <- model$specs$concomitant
    if (inherits(concom, "formula")) {
      X <- stats::model.matrix(concom, data = data)
    } else {
      f_str <- paste("~", paste(concom, collapse = " + "))
      X <- stats::model.matrix(stats::as.formula(f_str), data = data)
    }
  } else {
    X <- matrix(1, nrow = nrow(data), ncol = 1L)
  }
  coefs <- if (!is.null(model$concomitant_coefs)) model$concomitant_coefs else matrix(0, ncol(X), K - 1L)
  log_priors <- log(compute_priors(X, coefs) + 1e-300)

  cat_ld <- matrix(0, nrow = nrow(Y), ncol = K)
  if (has_cat) {
    for (k in seq_len(K))
      cat_ld[, k] <- eval_categorical_density(D, model$categorical_params[[k]])
  }

  # Observed marginal log-likelihood as function of continuous params
  obs_ll <- function(means_k, covs_k) {
    log_dens <- matrix(0, nrow = nrow(Y), ncol = K)
    for (k in seq_len(K)) {
      log_dens[, k] <- eval_continuous_density(Y, means_k[[k]], covs_k[[k]]) + cat_ld[, k]
    }
    log_joint <- log_priors + log_dens
    sum(apply(log_joint, 1, log_sum_exp))
  }

  means_list <- model$continuous_params$means
  covs_list  <- model$continuous_params$covariances
  ll_ct <- obs_ll(means_list, covs_list)

  mean_se <- list()
  cov_se  <- list()

  if (!requireNamespace("numDeriv", quietly = TRUE))
    stop("Package 'numDeriv' is required for continuous_se(). ",
         "Install via install.packages('numDeriv').")

  for (k in seq_len(K)) {
    mu_k    <- means_list[[k]]
    Sigma_k <- covs_list[[k]]

    # Class-specific mean SEs from full Hessian of observed-data log-likelihood
    ll_mu_k <- function(m) {
      means_temp       <- means_list
      means_temp[[k]]  <- m
      obs_ll(means_temp, covs_list)
    }
    H_mu  <- numDeriv::hessian(ll_mu_k, x = mu_k,
                               method.args = list(eps = 1e-4, d = 1e-4))
    V_mu  <- tryCatch(solve(-H_mu), error = function(e) matrix(NA_real_, d, d))
    se_mu <- suppressWarnings(sqrt(diag(V_mu)))
    se_mu[!is.finite(se_mu)] <- NA_real_
    names(se_mu) <- cont
    mean_se[[k]] <- se_mu

    # Diagonal-variance SEs from full Hessian over log-variance parameterization
    log_diag_k <- log(diag(Sigma_k))
    ll_logvar_k <- function(lv) {
      covs_temp        <- covs_list
      S_new            <- covs_temp[[k]]
      diag(S_new)      <- exp(lv)
      covs_temp[[k]]   <- S_new
      obs_ll(means_list, covs_temp)
    }
    H_lv  <- numDeriv::hessian(ll_logvar_k, x = log_diag_k,
                               method.args = list(eps = 1e-4, d = 1e-4))
    V_lv  <- tryCatch(solve(-H_lv), error = function(e) matrix(NA_real_, d, d))
    # Delta-rule: var(sigma2) = (sigma2)^2 * var(log sigma2)
    sigma2_diag <- diag(Sigma_k)
    se_diag <- suppressWarnings(sqrt(diag(V_lv)) * sigma2_diag)
    se_diag[!is.finite(se_diag)] <- NA_real_

    se_sig <- matrix(NA_real_, d, d)
    dimnames(se_sig) <- list(cont, cont)
    diag(se_sig) <- se_diag
    cov_se[[k]]  <- se_sig
  }

  list(mean_se = mean_se, cov_se = cov_se)
}
