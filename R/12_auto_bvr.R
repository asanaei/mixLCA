# ==============================================================================
# File: R/12_auto_bvr.R
# mixLCA - Automated model selection for number of classes and local dependence
# ==============================================================================

#' Automated Model Selection for mixLCA
#'
#' Three-phase automated model selection:
#' \describe{
#'   \item{Phase 1}{Select \emph{K} by BIC across \code{K_range} using
#'     diagonal covariance.}
#'   \item{Phase 2a}{Compare diagonal, penalized, and full covariance
#'     structures for continuous indicators at the selected \emph{K}.}
#'   \item{Phase 2b}{Vermunt-style BVR-guided specification search for
#'     categorical indicators: iteratively add direct effects for the
#'     pair with the largest bivariate residual, stopping when BIC no
#'     longer improves or \code{max_direct_effects} is reached.}
#' }
#'
#' Parallel execution of multiple starts is delegated to the user's active
#' \code{future::plan()} (set it in the session before calling this function).
#'
#' @param data Data frame.
#' @param continuous Character vector of continuous indicator names, or NULL.
#' @param categorical Character vector of categorical indicator names, or NULL.
#' @param concomitant Character vector or formula of concomitant predictors,
#'   or NULL.
#' @param K_range Integer vector of class counts to evaluate (e.g., 2:5).
#' @param max_direct_effects Integer: maximum number of categorical direct
#'   effects to add during Phase 2b (default 5).
#' @param bvr_threshold Numeric: minimum BVR chi-squared statistic to
#'   consider a direct effect (default 3.84, i.e. p < .05 for df = 1).
#' @param verbose Logical: print progress?
#' @param ... Additional arguments passed to \code{fit_lca} (e.g.,
#'   \code{spectral_rank}, \code{spectral_pool}).
#' @return An object of class \code{mixLCA} representing the final selected
#'   model, with an additional \code{auto_path} element recording the
#'   search trajectory.
#'
#' @examples
#' \dontrun{
#' set.seed(110)
#' data(voter_perceptions)
#' fit <- auto_bvr(
#'   data        = voter_perceptions,
#'   categorical = names(voter_perceptions),
#'   K_range     = 3,
#'   max_direct_effects = 4L,
#'   n_starts    = 3,
#'   verbose     = FALSE)
#' fit$auto_path$direct_effects
#' fit_indices(fit)$BIC
#' }
#' @export
auto_bvr <- function(data, continuous = NULL, categorical = NULL,
                     concomitant = NULL, K_range = 2:5,
                     max_direct_effects = 5L, bvr_threshold = 3.84,
                     verbose = TRUE, ...) {

  if (length(K_range) == 0) stop("K_range must not be empty.")

  dots <- list(...)
  user_n_starts <- if (!is.null(dots$n_starts)) dots$n_starts else 1L
  user_kmeans   <- if (!is.null(dots$kmeans_nstart)) dots$kmeans_nstart else 1L
  user_max_iter <- if (!is.null(dots$max_iter)) dots$max_iter else 500L
  user_tol      <- if (!is.null(dots$tol)) dots$tol else 1e-6
  dots$n_starts <- NULL
  dots$kmeans_nstart <- NULL
  dots$max_iter <- NULL
  dots$tol <- NULL
  dots$seed <- NULL
  dots$control <- NULL
  # Structural args injected by this wrapper into do.call(fit_lca, ...).
  # Clearing them prevents the "formal argument matched by multiple
  # actual arguments" fatal error if the user passes one of these
  # through `...`.
  dots$n_classes <- NULL
  dots$dependence <- NULL
  dots$penalty <- NULL
  dots$cat_direct_effects <- NULL
  dots$init_model <- NULL

  base_ctrl <- lca_control(
    max_iter = user_max_iter, tol = user_tol,
    n_starts = user_n_starts,
    kmeans_nstart = user_kmeans)
  warm_ctrl <- lca_control(
    max_iter = user_max_iter, tol = user_tol,
    n_starts = 1L,
    kmeans_nstart = user_kmeans)

  # ==== Phase 1: Select K ====
  if (verbose)
    message(sprintf("Phase 1: Selecting K across %d-%d",
                    min(K_range), max(K_range)))
  best_bic <- Inf
  best_K <- NA
  base_models <- list()

  for (i in seq_along(K_range)) {
    K <- K_range[i]
    concom_for_K <- if (K == 1) NULL else concomitant
    fit <- tryCatch({
      do.call(fit_lca, c(list(
        data = data, continuous = continuous, categorical = categorical,
        concomitant = concom_for_K, n_classes = K, dependence = "none",
        control = base_ctrl, verbose = FALSE), dots))
    }, error = function(e) NULL)

    if (is.null(fit) || !fit$convergence$converged) {
      if (verbose) message("  K = ", K, ": did not converge")
      next
    }

    base_models[[as.character(K)]] <- fit
    fi <- fit_indices(fit)
    if (verbose)
      message(sprintf("  K = %d: BIC = %.2f", K, fi$BIC))

    if (fi$BIC < best_bic) {
      best_bic <- fi$BIC
      best_K <- K
    }
  }

  if (is.na(best_K)) stop("No models converged in Phase 1.")
  if (verbose) message("  Selected K = ", best_K)

  current_model <- base_models[[as.character(best_K)]]
  current_bic   <- best_bic
  dep_list      <- list()

  # ==== Phase 2a: Continuous covariance structure ====
  if (!is.null(continuous) && length(continuous) >= 2) {
    if (verbose)
      message("Phase 2a: Evaluating continuous local dependence at K = ", best_K)

    pen_fit <- tryCatch({
      do.call(fit_lca, c(list(
        data = data, continuous = continuous, categorical = categorical,
        concomitant = concomitant, n_classes = best_K,
        dependence = "penalized", penalty = 0.1,
        control = warm_ctrl, verbose = FALSE,
        init_model = current_model), dots))
    }, error = function(e) NULL)

    if (!is.null(pen_fit) && pen_fit$convergence$converged) {
      pen_fi <- fit_indices(pen_fit)
      if (verbose) message(sprintf("  Penalized: BIC = %.2f", pen_fi$BIC))
      if (pen_fi$BIC < current_bic) {
        current_model <- pen_fit
        current_bic   <- pen_fi$BIC
      }
    }

    full_fit <- tryCatch({
      do.call(fit_lca, c(list(
        data = data, continuous = continuous, categorical = categorical,
        concomitant = concomitant, n_classes = best_K,
        dependence = "full",
        control = warm_ctrl, verbose = FALSE,
        init_model = current_model), dots))
    }, error = function(e) NULL)

    if (!is.null(full_fit) && full_fit$convergence$converged) {
      full_fi <- fit_indices(full_fit)
      if (verbose) message(sprintf("  Full: BIC = %.2f", full_fi$BIC))
      if (full_fi$BIC < current_bic) {
        current_model <- full_fit
        current_bic   <- full_fi$BIC
      }
    }
  }

  # ==== Phase 2b: Categorical BVR-guided specification search ====
  has_spectral <- is.numeric(dots$spectral_rank) && any(dots$spectral_rank > 0L)
  if (!is.null(categorical) && length(categorical) >= 2 && !has_spectral) {
    if (verbose)
      message("Phase 2b: Categorical BVR specification search at K = ", best_K)

    dep_dependence <- current_model$specs$dependence

    for (step in seq_len(max_direct_effects)) {
      bvr_tab <- tryCatch(
        bvr_categorical(current_model, data),
        error = function(e) NULL)
      if (is.null(bvr_tab) || nrow(bvr_tab) == 0) break

      existing <- vapply(dep_list, paste, "", collapse = "~")

      candidate <- NULL
      for (r in seq_len(nrow(bvr_tab))) {
        pair_key  <- paste(bvr_tab$var1[r], bvr_tab$var2[r], sep = "~")
        pair_key2 <- paste(bvr_tab$var2[r], bvr_tab$var1[r], sep = "~")
        if (bvr_tab$bvr[r] < bvr_threshold) break
        if (!(pair_key %in% existing) && !(pair_key2 %in% existing)) {
          candidate <- c(bvr_tab$var1[r], bvr_tab$var2[r])
          candidate_bvr <- bvr_tab$bvr[r]
          break
        }
      }
      if (is.null(candidate)) break

      trial_deps <- c(dep_list, list(candidate))
      trial_fit <- tryCatch({
        do.call(fit_lca, c(list(
          data = data, continuous = continuous, categorical = categorical,
          concomitant = concomitant, n_classes = best_K,
          dependence = dep_dependence,
          control = warm_ctrl, verbose = FALSE,
          cat_direct_effects = trial_deps,
          init_model = current_model), dots))
      }, error = function(e) NULL)

      if (is.null(trial_fit) || !trial_fit$convergence$converged) {
        if (verbose) message("  ", candidate[1], " -> ", candidate[2],
                             ": did not converge, skipping")
        next
      }

      trial_fi <- fit_indices(trial_fit)
      if (verbose)
        message(sprintf("  %s -> %s: BVR = %.1f, BIC = %.2f",
                        candidate[1], candidate[2],
                        candidate_bvr, trial_fi$BIC))

      if (trial_fi$BIC < current_bic) {
        dep_list      <- trial_deps
        current_model <- trial_fit
        current_bic   <- trial_fi$BIC
      } else {
        if (verbose) message("  BIC did not improve; stopping search.")
        break
      }
    }
  }

  if (verbose)
    message("  Final BIC = ", round(current_bic, 2),
            " (dependence: ", current_model$specs$dependence,
            if (length(dep_list) > 0)
              paste0(", ", length(dep_list), " direct effect(s)"),
            ")")

  current_model$auto_path <- list(
    best_K            = best_K,
    final_bic         = current_bic,
    direct_effects    = dep_list
  )
  return(current_model)
}
