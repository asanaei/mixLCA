# ==============================================================================
# File: R/06_diagnostics.R
# mixLCA - model fit indices, bivariate residual analysis,
# classification diagnostics, multi-model comparison.
# ==============================================================================

#' Compute Fit Indices for a mixLCA Model
#'
#' Returns AIC, BIC, sample-adjusted BIC, relative entropy, and the
#' Integrated Classification Likelihood criterion (ICL).
#'
#' @param model A \code{mixLCA} object.
#' @return Named list of fit statistics.
#'
#' @examples
#' \donttest{
#' data(voter_perceptions)
#' fit <- fit_lca(voter_perceptions,
#'                categorical = names(voter_perceptions),
#'                n_classes   = 3,
#'                control     = lca_control(n_starts = 2, seed = 110),
#'                verbose     = FALSE)
#' fi <- fit_indices(fit)
#' fi[c("log_lik", "n_params", "AIC", "BIC", "aBIC", "entropy", "ICL")]
#' }
#' @export
fit_indices <- function(model) {
  if (!inherits(model, "mixLCA"))
    stop("`model` must be a mixLCA object.")

  ll <- model$log_lik
  np <- model$n_params
  N  <- if (!is.null(model$n_obs_effective) && !is.na(model$n_obs_effective))
          model$n_obs_effective else model$n_obs
  K  <- model$n_classes

  aic  <- -2 * ll + 2 * np
  bic  <- -2 * ll + np * log(N)
  abic <- -2 * ll + np * log((N + 2) / 24)

  # Relative entropy (0 = no separation, 1 = crisp assignment)
  post      <- model$posteriors
  post_safe <- pmax(post, 1e-15)
  ent_raw   <- -sum(post_safe * log(post_safe))
  entropy   <- if (K > 1L) 1 - ent_raw / (N * log(K)) else 1

  # ICL = BIC penalised by classification uncertainty
  icl <- bic + 2 * ent_raw

  # Composite-likelihood flag
  is_composite <- isTRUE(any(model$specs$spectral_rank > 0L))

  out <- list(
    log_lik      = ll,
    n_params     = np,
    AIC          = aic,
    BIC          = bic,
    aBIC         = abic,
    entropy      = entropy,
    ICL          = icl,
    is_composite = is_composite
  )
  if (is_composite) {
    out$BIC_note <- "naive (composite likelihood); prefer cross-validation for rank selection"
  }

  out
}

#' Bivariate Residual Covariance Matrix
#'
#' Computes the residual covariance between continuous indicators after
#' subtracting the model-implied mixture covariance. Large residuals
#' signal local dependence violations not captured by the current
#' model. Use \code{plot_bvr()} to visualise the residual network.
#'
#' @param model A \code{mixLCA} object.
#' @param data Data frame used for estimation.
#' @return Named matrix of residual covariances among continuous
#'   indicators.
#'
#' @examples
#' \donttest{
#' data(health_screening)
#' fit <- fit_lca(health_screening,
#'                continuous = c("marker_1","marker_2","marker_3","marker_4"),
#'                n_classes  = 2,
#'                control    = lca_control(n_starts = 2, seed = 110),
#'                verbose    = FALSE)
#' round(bvr(fit, health_screening), 3)
#' }
#' @export
bvr <- function(model, data) {
  if (!inherits(model, "mixLCA"))
    stop("`model` must be a mixLCA object.")

  cont <- model$specs$continuous
  if (is.null(cont))
    stop("bvr() requires continuous indicators.")

  Y <- as.matrix(data[, cont, drop = FALSE])
  d <- ncol(Y)
  K <- model$n_classes

  S_obs    <- stats::cov(Y, use = "pairwise.complete.obs")
  pi_k     <- colMeans(model$posteriors)
  mu_global <- numeric(d)
  for (k in seq_len(K))
    mu_global <- mu_global +
      pi_k[k] * model$continuous_params$means[[k]]

  S_imp <- matrix(0, d, d)
  for (k in seq_len(K)) {
    delta <- model$continuous_params$means[[k]] - mu_global
    S_imp <- S_imp + pi_k[k] *
      (model$continuous_params$covariances[[k]] + tcrossprod(delta))
  }

  resid <- S_obs - S_imp
  rownames(resid) <- cont
  colnames(resid) <- cont
  resid
}

#' Bivariate Residual Significance Tests
#'
#' Converts residual covariances to approximate chi-squared statistics
#' (df = 1) for pairwise local dependence testing.
#'
#' @param model A \code{mixLCA} object.
#' @param data Data frame used for estimation.
#' @return Data frame with columns: \code{var1}, \code{var2},
#'   \code{residual_cov}, \code{chi_sq}, \code{p_value}.
#'
#' @examples
#' \donttest{
#' data(health_screening)
#' fit <- fit_lca(health_screening,
#'                continuous = c("marker_1","marker_2","marker_3","marker_4"),
#'                n_classes  = 2,
#'                control    = lca_control(n_starts = 2, seed = 110),
#'                verbose    = FALSE)
#' bvr_tests(fit, health_screening)
#' }
#' @export
bvr_tests <- function(model, data) {
  bvr_mat  <- bvr(model, data)
  cont <- model$specs$continuous
  N    <- model$n_obs
  d    <- length(cont)

  rows <- list()
  idx  <- 1L
  for (i in seq_len(d - 1L)) {
    for (j in (i + 1L):d) {
      r_ij  <- bvr_mat[i, j]
      s_i   <- sqrt(stats::var(data[[cont[i]]], na.rm = TRUE))
      s_j   <- sqrt(stats::var(data[[cont[j]]], na.rm = TRUE))
      r_std <- r_ij / (s_i * s_j + 1e-15)
      chi2  <- N * r_std^2
      pval  <- stats::pchisq(chi2, df = 1, lower.tail = FALSE)
      rows[[idx]] <- data.frame(
        var1         = cont[i],
        var2         = cont[j],
        residual_cov = r_ij,
        chi_sq       = chi2,
        p_value      = pval,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  do.call(rbind, rows)
}

#' Classification Table for a mixLCA Model
#'
#' Cross-tabulates modal class assignments with average posterior
#' probabilities. Rows = modal assignment, columns = average P(class k).
#' Each row sums to 1. Diagonal dominance indicates good class separation.
#'
#' @param model A \code{mixLCA} object.
#' @return K x K matrix: rows = modal assignment, columns = average
#'   P(class k).
#'
#' @examples
#' \donttest{
#' data(voter_perceptions)
#' fit <- fit_lca(voter_perceptions,
#'                categorical = names(voter_perceptions),
#'                n_classes   = 3,
#'                control     = lca_control(n_starts = 2, seed = 110),
#'                verbose     = FALSE)
#' round(class_table(fit), 3)
#' }
#' @export
class_table <- function(model) {
  if (!inherits(model, "mixLCA"))
    stop("`model` must be a mixLCA object.")

  posteriors <- model$posteriors
  K <- ncol(posteriors)
  modal <- apply(posteriors, 1, which.max)
  E <- matrix(0, K, K)

  for (w in seq_len(K)) {
    idx <- which(modal == w)
    if (length(idx) > 0L) {
      E[w, ] <- colMeans(posteriors[idx, , drop = FALSE])
    } else {
      E[w, w] <- 1
    }
  }
  E
}

#' Compare Multiple mixLCA Models
#'
#' Assembles a summary table of fit indices across a list of fitted
#' models, ordered by number of classes.
#'
#' @param ... One or more \code{mixLCA} objects, or a single named list
#'   of them.
#' @return Data frame of comparative fit statistics.
#'
#' @examples
#' \donttest{
#' data(voter_perceptions)
#' cat_items <- names(voter_perceptions)
#' fits <- lapply(2:4, function(K)
#'   fit_lca(voter_perceptions, categorical = cat_items, n_classes = K,
#'           control = lca_control(n_starts = 2, seed = 110),
#'           verbose = FALSE))
#' names(fits) <- paste0("K", 2:4)
#' compare_models(fits)
#' }
#' @export
compare_models <- function(...) {
  models <- list(...)
  if (length(models) == 1L && is.list(models[[1L]]) &&
      !inherits(models[[1L]], "mixLCA"))
    models <- models[[1L]]

  spectral_maxes <- sapply(models, function(m) {
    sr <- m$specs$spectral_rank
    if (is.null(sr)) 0L else max(sr)
  })
  if (length(unique(spectral_maxes)) > 1L) {
    warning("Comparing models with different spectral ranks (d) using standard AIC/BIC is heuristic. ",
            "Standard criteria typically under-penalize composite likelihoods (d > 0). ",
            "Use out-of-sample cross-validation to formally select rank.", call. = FALSE)
  }

  rows <- lapply(models, function(m) {
    fi <- fit_indices(m)
    data.frame(
      K        = m$n_classes,
      LL       = fi$log_lik,
      n_params = fi$n_params,
      AIC      = fi$AIC,
      BIC      = fi$BIC,
      aBIC     = fi$aBIC,
      entropy  = fi$entropy,
      ICL      = fi$ICL,
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  tab[order(tab$K), ]
}
