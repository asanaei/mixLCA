# ==============================================================================
# File: R/15_auto_spectral.R
# mixLCA - Adaptive greedy search for class-specific spectral ranks
# ==============================================================================

#' Extract next unmodeled eigenvalue per class
#'
#' For each class with rank below the ceiling, retrieves the next
#' eigenvalue from the conditional Burt matrix spectrum. If the model
#' already stores the full spectrum (as it does after at least one EM
#' cycle with SLD active), the value is read directly; otherwise, the
#' weighted residual Burt matrix is computed from the posteriors.
#' @keywords internal
get_unmodeled_eigenvalues <- function(model, data, ranks, max_rank_per_class) {
  K <- model$n_classes
  eigs <- rep(-Inf, K)
  D <- data[, model$specs$categorical, drop = FALSE]
  spec_data <- prep_spectral_data(D)
  post <- model$posteriors

  for (k in seq_len(K)) {
    if (ranks[k] >= max_rank_per_class) next

    # Fast path: read stored spectrum, but only trust it when this class
    # has rank > 0 (so the stored eigenvalues come from a real decomposition).
    # Classes at rank 0 carry a placeholder zero-vector, which would falsely
    # tell us "no remaining eigenvalue" and stop the search.
    if (ranks[k] > 0L &&
        !is.null(model$cat_spectral_params) &&
        length(model$cat_spectral_params$spectra) >= k &&
        length(model$cat_spectral_params$spectra[[k]]) > ranks[k]) {
      val <- model$cat_spectral_params$spectra[[k]][ranks[k] + 1L]
      if (is.numeric(val) && !is.na(val)) { eigs[k] <- val; next }
    }

    # Slow path: compute from posteriors
    wk <- post[, k]
    W_sum <- sum(wk)
    if (W_sum > 1e-15) {
      pi_c <- flatten_cat_params(model$categorical_params[[k]],
                                 spec_data$item_indices, spec_data$C)
      E_c <- matrix(pi_c, nrow = nrow(D), ncol = spec_data$C, byrow = TRUE)
      Z_imp <- spec_data$Z
      Z_imp[spec_data$Z_mis] <- E_c[spec_data$Z_mis]
      R_c <- Z_imp - E_c
      R_w <- sweep(R_c, 1, sqrt(wk), "*")
      Sigma_k <- crossprod(R_w) / W_sum
      Sigma_k <- Sigma_k * spec_data$M
      eig <- eigen((Sigma_k + t(Sigma_k)) / 2, symmetric = TRUE, only.values = TRUE)
      val <- eig$values[ranks[k] + 1L]
      if (is.numeric(val) && !is.na(val)) eigs[k] <- val
    }
  }
  eigs
}

#' Automated Class-Specific Spectral Rank Selection
#'
#' Evaluates the necessity of Spectral Local Dependence (SLD) iteratively.
#' Begins with a base model (all ranks zero or a user-supplied fitted
#' model) and conducts a forward-stepwise greedy search. At each step
#' the unmodeled eigenvalue spectrum is inspected per class, and the
#' class exhibiting the largest residual eigenvalue is targeted for a
#' rank increment. If the resulting model improves the information
#' criterion, the increment is accepted; otherwise the search terminates.
#'
#' Because each candidate is hot-started from the current model, the
#' latent class definitions remain anchored, preventing latent-class
#' drift. Classes that do not exhibit residual dependence retain
#' \code{d = 0}, avoiding unnecessary parameter proliferation.
#'
#' @param data Data frame.
#' @param continuous Character vector or NULL.
#' @param categorical Character vector (at least 2 required).
#' @param concomitant Character vector or NULL.
#' @param n_classes Integer: number of latent classes.
#' @param max_rank_per_class Integer: max spectral rank for any class.
#' @param max_total_rank Integer: max sum of ranks across all classes.
#' @param criterion Character: information criterion to minimize
#'   (\code{"BIC"}, \code{"aBIC"}, \code{"AIC"}, or \code{"ICL"}).
#' @param base_model Optional pre-fitted \code{mixLCA} object. If NULL,
#'   an independence model (d = 0) is fitted automatically.
#' @param seed Integer random seed for the base model.
#' @param verbose Logical: print search trajectory.
#' @param ... Additional arguments forwarded to \code{fit_lca}
#'   (e.g., \code{max_iter}, \code{tol}, \code{dependence}).
#'
#' @return A \code{mixLCA} object of the selected configuration, with
#'   an additional element \code{auto_spectral_path} documenting the
#'   search history (data frame with columns \code{step},
#'   \code{class_incremented}, \code{ranks}, and the criterion value).
#'
#' @examples
#' \dontrun{
#' data(voter_perceptions)
#' fit <- auto_sld(
#'   data        = voter_perceptions,
#'   categorical = names(voter_perceptions),
#'   n_classes   = 3,
#'   max_rank_per_class = 3L,
#'   criterion   = "BIC",
#'   seed        = 110L,
#'   verbose     = FALSE)
#' fit$specs$spectral_rank      # class-specific ranks selected
#' fit$auto_spectral_path        # accepted increments
#' fit_indices(fit)$BIC
#' }
#' @export
auto_sld <- function(data, continuous = NULL, categorical = NULL,
                     concomitant = NULL, n_classes = 2L,
                     max_rank_per_class = 3L, max_total_rank = 10L,
                     criterion = c("BIC", "aBIC", "AIC", "ICL"),
                     base_model = NULL, seed = 110L,
                     verbose = TRUE, ...) {

  if (is.null(categorical) || length(categorical) < 2L)
    stop("Spectral Local Dependence requires at least 2 categorical indicators.")

  criterion <- match.arg(criterion)
  K <- as.integer(n_classes)

  dots <- list(...)
  user_max_iter <- if (!is.null(dots$max_iter)) dots$max_iter else 500L
  user_tol      <- if (!is.null(dots$tol)) dots$tol else 1e-6
  user_kmeans   <- if (!is.null(dots$kmeans_nstart)) dots$kmeans_nstart else 1L
  dots$n_starts <- NULL
  dots$max_iter <- NULL
  dots$tol <- NULL
  dots$seed <- NULL
  dots$kmeans_nstart <- NULL
  dots$control <- NULL

  ctrl <- lca_control(
    max_iter = user_max_iter, tol = user_tol,
    n_starts = 1L, seed = seed,
    kmeans_nstart = user_kmeans)

  get_ic <- function(model, crit) fit_indices(model)[[crit]]

  if (verbose) {
    message("======================================================")
    message("Adaptive Spectral Local Dependence Search (Auto-SLD)")
    message("======================================================")
  }

  # Fit or validate the baseline model
  if (is.null(base_model)) {
    if (verbose) message("Fitting baseline model (d = 0 for all classes)...")
    current_model <- tryCatch(
      do.call(fit_lca, c(list(
        data = data, continuous = continuous,
        categorical = categorical, concomitant = concomitant,
        n_classes = K, spectral_rank = rep(0L, K),
        spectral_pool = FALSE, control = ctrl,
        verbose = FALSE), dots)),
      error = function(e) stop("Baseline model failed: ", e$message))
  } else {
    if (!inherits(base_model, "mixLCA"))
      stop("base_model must be a mixLCA object.")
    if (base_model$n_classes != K)
      stop("base_model has different number of classes than n_classes.")
    current_model <- base_model
  }

  if (!current_model$convergence$converged)
    stop("Base model failed to converge.")

  current_ranks <- current_model$specs$spectral_rank
  if (is.null(current_ranks)) current_ranks <- rep(0L, K)
  if (length(current_ranks) == 1L) current_ranks <- rep(current_ranks, K)

  current_ic <- get_ic(current_model, criterion)
  if (verbose) message(sprintf("Baseline %s: %.2f\n", criterion, current_ic))

  search_path <- list(
    data.frame(step = 0L, class_incremented = NA_integer_,
               ranks = paste(current_ranks, collapse = ", "),
               ic = current_ic, stringsAsFactors = FALSE))
  step <- 1L

  while (sum(current_ranks) < max_total_rank) {
    if (verbose) message(sprintf("--- Step %d ---", step))

    # O(1) targeting: identify the class with the largest unmodeled eigenvalue
    unmodeled_eigs <- get_unmodeled_eigenvalues(
      current_model, data, current_ranks, max_rank_per_class)
    best_candidate_k <- as.integer(which.max(unmodeled_eigs))
    max_next_eig <- unmodeled_eigs[best_candidate_k]

    if (length(best_candidate_k) == 0L || is.na(best_candidate_k) ||
        max_next_eig < 1e-5) {
      if (verbose)
        message("=> STOP: No valid remaining eigenvalues to model. Search complete.\n")
      break
    }

    trial_ranks <- current_ranks
    trial_ranks[best_candidate_k] <- trial_ranks[best_candidate_k] + 1L

    if (verbose)
      message(sprintf("  -> Targeting Class %d (next unmodeled eigenvalue = %.4f)",
                      best_candidate_k, max_next_eig))

    trial_model <- tryCatch({
      do.call(fit_lca, c(list(
        data = data, continuous = continuous,
        categorical = categorical, concomitant = concomitant,
        n_classes = K, spectral_rank = trial_ranks,
        spectral_pool = FALSE, control = ctrl,
        verbose = FALSE, init_model = current_model), dots))
    }, error = function(e) NULL)

    # Treat a candidate as acceptable if either the converged flag is set
    # or the log-likelihood has essentially flattened (mean of the last
    # three iterations' increments below 1e-3, which is loose relative to
    # BIC differences but tight relative to the iteration cap).
    soft_converged <- function(m) {
      if (is.null(m)) return(FALSE)
      if (isTRUE(m$convergence$converged)) return(TRUE)
      hist <- m$convergence$ll_history
      if (length(hist) < 4L) return(FALSE)
      mean(abs(diff(utils::tail(hist, 4L)))) < 1e-3
    }

    if (!is.null(trial_model) && soft_converged(trial_model)) {
      trial_ic <- get_ic(trial_model, criterion)
      if (verbose) {
        tag <- if (trial_ic < current_ic) " (improved)" else ""
        message(sprintf("  Test ranks [%s]: %s = %.2f%s",
                        paste(trial_ranks, collapse = ", "),
                        criterion, trial_ic, tag))
      }
      if (trial_ic < current_ic) {
        current_ranks <- trial_ranks
        current_model <- trial_model
        current_ic    <- trial_ic
        if (verbose)
          message(sprintf("=> ACCEPTED: rank +1 for Class %d. New %s: %.2f\n",
                          best_candidate_k, criterion, current_ic))

        search_path[[step + 1L]] <- data.frame(
          step = step, class_incremented = best_candidate_k,
          ranks = paste(current_ranks, collapse = ", "),
          ic = current_ic, stringsAsFactors = FALSE)
        step <- step + 1L
      } else {
        if (verbose)
          message("=> STOP: Information criterion did not improve. Search complete.\n")
        break
      }
    } else {
      if (verbose)
        message("=> STOP: Candidate model failed to converge. Search complete.\n")
      break
    }
  }

  if (verbose) {
    message("======================================================")
    message(sprintf("Final ranks : [%s]", paste(current_ranks, collapse = ", ")))
    message(sprintf("Final %s    : %.2f", criterion, current_ic))
  }

  current_model$auto_spectral_path <- do.call(rbind, search_path)
  current_model
}
