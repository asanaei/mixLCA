# ==============================================================================
# File: R/17_tidyclust_bridge.R
# mixLCA - bridge functions between tidyclust interface and fit_lca/predict,
# plus broom generics (tidy, glance) for tidymodels integration.
# ==============================================================================

#' Fit Bridge for tidyclust
#'
#' Translates the tidyclust \code{fit()} call into a
#' \code{\link{fit_lca}()} call.  Column roles (continuous vs.
#' categorical) are specified through engine arguments; when omitted
#' they are inferred from column types.  Concomitant variables are
#' excluded from auto-detection to prevent double-counting.
#'
#' @param x Data frame of predictors (intercept already removed by
#'   tidyclust encoding).
#' @param num_clusters Integer number of latent classes.
#' @param dependence Covariance structure passed to \code{fit_lca()}.
#' @param ... Additional engine arguments forwarded to
#'   \code{fit_lca()} (e.g. \code{categorical}, \code{continuous},
#'   \code{concomitant}, \code{control}).
#' @return A fitted \code{mixLCA} object.
#' @keywords internal
#' @export
.lca_clust_fit_mixLCA <- function(x, num_clusters, dependence = "full",
                                  penalty = NULL, spectral_rank = NULL,
                                  ...) {
  dots <- list(...)

  continuous  <- dots$continuous
  categorical <- dots$categorical
  concomitant <- dots$concomitant
  dots$continuous  <- NULL
  dots$categorical <- NULL

  if (is.null(continuous) && is.null(categorical)) {
    concom_vars <- if (inherits(concomitant, "formula")) {
      all.vars(concomitant)
    } else {
      concomitant
    }

    ind_vars <- setdiff(names(x), concom_vars)
    if (length(ind_vars) == 0L)
      stop("No manifest indicators remain after excluding concomitant variables.")

    is_num      <- vapply(x[ind_vars], is.numeric, logical(1L))
    continuous  <- ind_vars[is_num]
    categorical <- ind_vars[!is_num]
    if (length(continuous)  == 0L) continuous  <- NULL
    if (length(categorical) == 0L) categorical <- NULL
  }

  fit_args <- c(
    list(
      data          = x,
      continuous    = continuous,
      categorical   = categorical,
      n_classes     = as.integer(num_clusters),
      dependence    = dependence,
      penalty       = penalty,
      spectral_rank = spectral_rank,
      verbose       = FALSE
    ),
    dots
  )

  fit_args <- fit_args[!vapply(fit_args, is.null, logical(1L))]

  do.call(fit_lca, fit_args)
}

#' Predict Bridge for tidyclust (cluster assignments)
#'
#' Returns a factor of cluster labels. tidyclust's
#' \code{predict.cluster_fit} wraps this into the standard
#' \code{.pred_cluster} tibble column.
#'
#' @param object A fitted \code{mixLCA} object.
#' @param new_data Data frame of new observations.
#' @return Factor of cluster assignments.
#' @keywords internal
#' @export
.lca_clust_predict_mixLCA <- function(object, new_data) {
  cls <- predict.mixLCA(object, newdata = new_data, type = "class")
  K <- object$n_classes
  factor(cls, levels = seq_len(K),
         labels = paste0("Cluster_", seq_len(K)))
}

#' Predict Bridge for tidyclust (class probabilities)
#'
#' Returns a data frame of posterior class probabilities with
#' columns named \code{.pred_Cluster_1}, \code{.pred_Cluster_2},
#' etc., as required by the tidyclust predict interface.
#'
#' @param object A fitted \code{mixLCA} object.
#' @param new_data Data frame of new observations.
#' @return Data frame of posterior probabilities.
#' @keywords internal
#' @export
.lca_clust_predict_prob_mixLCA <- function(object, new_data) {
  probs <- predict.mixLCA(object, newdata = new_data, type = "prob")
  K <- object$n_classes
  colnames(probs) <- paste0(".pred_Cluster_", seq_len(K))
  as.data.frame(probs)
}

# ------------------------------------------------------------------------------
# broom generics: tidy and glance
# ------------------------------------------------------------------------------

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' Tidy a mixLCA Model
#'
#' Summarizes each latent class as one row, with size and mixing
#' proportion.
#'
#' @param x A \code{mixLCA} object.
#' @param ... Unused.
#' @return A data frame with columns \code{cluster}, \code{size},
#'   \code{proportion}.
#'
#' @examples
#' \donttest{
#' data(voter_perceptions)
#' fit <- fit_lca(voter_perceptions,
#'                categorical = names(voter_perceptions),
#'                n_classes   = 3,
#'                control     = lca_control(n_starts = 2),
#'                verbose     = FALSE)
#' tidy(fit)
#' }
#' @export
tidy.mixLCA <- function(x, ...) {
  props <- colMeans(x$posteriors)
  data.frame(
    cluster    = paste0("Cluster_", seq_len(x$n_classes)),
    size       = as.integer(round(colSums(x$posteriors))),
    proportion = props,
    stringsAsFactors = FALSE
  )
}

#' Glance at a mixLCA Model
#'
#' Returns a one-row data frame of fit statistics suitable for
#' model comparison tables and \code{tune_cluster()} logging.
#'
#' @param x A \code{mixLCA} object.
#' @param ... Unused.
#' @return A one-row data frame with columns \code{n_classes},
#'   \code{logLik}, \code{AIC}, \code{BIC}, \code{aBIC},
#'   \code{entropy}, \code{ICL}, \code{nobs}, \code{n_params},
#'   \code{converged}, \code{iterations}.
#'
#' @examples
#' \donttest{
#' data(voter_perceptions)
#' fit <- fit_lca(voter_perceptions,
#'                categorical = names(voter_perceptions),
#'                n_classes   = 3,
#'                control     = lca_control(n_starts = 2),
#'                verbose     = FALSE)
#' glance(fit)
#' }
#' @export
glance.mixLCA <- function(x, ...) {
  fi <- fit_indices(x)
  data.frame(
    n_classes  = x$n_classes,
    logLik     = fi$log_lik,
    AIC        = fi$AIC,
    BIC        = fi$BIC,
    aBIC       = fi$aBIC,
    entropy    = fi$entropy,
    ICL        = fi$ICL,
    nobs       = nobs(x),
    n_params   = fi$n_params,
    converged  = x$convergence$converged,
    iterations = x$convergence$iterations,
    stringsAsFactors = FALSE
  )
}
