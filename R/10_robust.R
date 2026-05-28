# ==============================================================================
# File: R/10_robust.R
# mixLCA - multi-start estimation and model enumeration.
# ==============================================================================

#' Multi-Start EM Estimation for mixLCA
#'
#' Runs the EM algorithm from multiple random starting points (each
#' seeded deterministically from \code{base_seed + s}) and returns the
#' solution with the highest terminal log-likelihood. Parallel execution
#' is delegated to the user's active \code{future::plan()}; this function
#' never alters the global plan or the global random seed.
#'
#' @param data Data frame.
#' @param continuous Character vector or NULL.
#' @param categorical Character vector or NULL.
#' @param concomitant Character vector or formula or NULL.
#' @param n_classes Integer >= 2.
#' @param dependence Character.
#' @param penalty Numeric.
#' @param max_iter Integer.
#' @param tol Numeric.
#' @param n_starts Number of random starting configurations.
#' @param base_seed Base random seed.
#' @param verbose Logical: print per-start progress?
#' @param cat_direct_effects List of direct effect pairs, or NULL.
#' @param spectral_rank Integer: SLD rank (0 = disabled).
#' @param spectral_pool Logical: pool Burt matrices across classes.
#' @param init_model Optional prior \code{mixLCA} object for warm-start
#'   (used for the first start only).
#' @param kmeans_nstart Integer: random starts for internal k-means.
#' @return The best-fitting \code{mixLCA} object.
#' @keywords internal
run_em_robust <- function(data, continuous = NULL, categorical = NULL,
                          concomitant = NULL, n_classes = 2L,
                          dependence = "full", penalty = 0,
                          max_iter = 500L, tol = 1e-6,
                          n_starts = 1L, base_seed = 110L,
                          verbose = TRUE,
                          cat_direct_effects = NULL,
                          spectral_rank = 0L, spectral_pool = FALSE,
                          init_model = NULL, kmeans_nstart = 1L) {

  best_model <- NULL
  best_ll    <- -Inf

  fit_start <- function(s) {
    seed_s <- base_seed + s
    candidate <- tryCatch(
      run_em(
        data        = data,
        continuous  = continuous,
        categorical = categorical,
        concomitant = concomitant,
        n_classes   = n_classes,
        dependence  = dependence,
        penalty     = penalty,
        max_iter    = max_iter,
        tol         = tol,
        seed        = seed_s,
        cat_direct_effects = cat_direct_effects,
        spectral_rank = spectral_rank,
        spectral_pool = spectral_pool,
        init_model  = if (s == 1L) init_model else NULL,
        kmeans_nstart = kmeans_nstart
      ),
      error = function(e) {
        if (verbose)
          message("  Start ", s, " failed: ", conditionMessage(e))
        NULL
      }
    )
    return(candidate)
  }

  if (n_starts > 1L && requireNamespace("future.apply", quietly = TRUE)) {
    # Delegate to the user's active future::plan(); future.seed = TRUE
    # safely manages parallel RNG without modifying .Random.seed.
    results <- future.apply::future_lapply(seq_len(n_starts), fit_start,
                                           future.seed = TRUE)
  } else {
    results <- lapply(seq_len(n_starts), fit_start)
  }

  for (s in seq_len(n_starts)) {
    candidate <- results[[s]]

    if (!is.null(candidate)) {
      if (verbose)
        message("  Start ", s, ": LL = ",
                round(candidate$log_lik, 4),
                if (candidate$convergence$converged) " (converged)"
                else paste0(" (", candidate$convergence$iterations,
                            " iters)"))

      if (candidate$log_lik > best_ll) {
        best_ll    <- candidate$log_lik
        best_model <- candidate
      }
    }
  }

  if (is.null(best_model))
    stop("All ", n_starts, " starts failed. ",
         "Examine the data for severe collinearity or rank deficiency.")

  if (verbose)
    message("  Best LL across starts: ", round(best_ll, 4))

  best_model
}

#' Enumerate mixLCA Models Across Class Counts
#'
#' Fits \code{fit_lca()} for each value in \code{k_range} and returns
#' all fitted objects together with a model comparison table. To run
#' the starts in parallel, set a \code{future::plan()} in your session
#' before calling this function.
#'
#' @param data Data frame.
#' @param continuous Character vector or NULL.
#' @param categorical Character vector or NULL.
#' @param concomitant Character vector or formula or NULL.
#' @param k_range Integer vector of class counts to estimate, e.g.\cr
#'   \code{2:5}.
#' @param dependence Character.
#' @param penalty Numeric.
#' @param n_starts Integer.
#' @param max_iter Integer.
#' @param tol Numeric.
#' @param spectral_rank Integer: SLD rank (0 = disabled).
#' @param spectral_pool Logical: pool Burt matrices across classes.
#' @param verbose Logical.
#' @param kmeans_nstart Integer: random starts for internal k-means (default 1).
#' @return List with:
#'   \describe{
#'     \item{\code{$models}}{Named list of \code{mixLCA} objects.}
#'     \item{\code{$comparison}}{Data frame from \code{compare_models()}.}
#'   }
#' @export
enumerate_lca <- function(data, continuous = NULL, categorical = NULL,
                          concomitant = NULL, k_range = 2:4,
                          dependence = "full", penalty = 0,
                          n_starts = 1L, max_iter = 500L, tol = 1e-6,
                          spectral_rank = 0L, spectral_pool = FALSE,
                          verbose = TRUE, kmeans_nstart = 1L) {

  validate_inputs(data, continuous, categorical, concomitant,
                  n_classes = min(k_range), dependence, penalty)

  models <- list()
  for (K in k_range) {
    if (verbose) message("Estimating K = ", K, " ...")
    models[[paste0("K", K)]] <- run_em_robust(
      data        = data,
      continuous  = continuous,
      categorical = categorical,
      concomitant = concomitant,
      n_classes   = K,
      dependence  = dependence,
      penalty     = penalty,
      max_iter    = max_iter,
      tol         = tol,
      n_starts    = n_starts,
      base_seed   = 110L,
      spectral_rank = spectral_rank,
      spectral_pool = spectral_pool,
      verbose     = verbose,
      kmeans_nstart = kmeans_nstart
    )
  }

  tab <- compare_models(models)
  if (verbose) {
    cat("\nModel Comparison (mixLCA):\n")
    print(tab, row.names = FALSE)
    cat("\n")
  }

  list(models = models, comparison = tab)
}
