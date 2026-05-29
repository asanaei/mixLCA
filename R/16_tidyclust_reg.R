# ==============================================================================
# File: R/16_tidyclust_reg.R
# mixLCA - tidyclust model specification, engine registration, and
# .onLoad hook.
# ==============================================================================

#' Latent Class Clustering
#'
#' @description
#'
#' \code{lca_clust()} defines a latent class model for clustering via
#' finite mixture estimation.
#'
#' This specification is designed for use with the
#' \href{https://tidyclust.tidymodels.org}{tidyclust} framework.
#' The engine-specific details for the \pkg{mixLCA} engine are
#' documented below.
#'
#' @param mode A single character string for the type of model. The
#'   only possible value for this model is \code{"partition"}.
#' @param engine A single character string specifying the
#'   computational engine. Currently only \code{"mixLCA"}.
#' @param num_clusters Positive integer, the number of latent classes
#'   (required).
#' @param dependence Character controlling the covariance structure of
#'   continuous indicators: \code{"none"} (diagonal),
#'   \code{"full"} (unrestricted), or \code{"penalized"}
#'   (graphical-lasso).
#' @param penalty Non-negative numeric penalty for graphical-lasso
#'   estimation (used when \code{dependence = "penalized"}), or
#'   \code{"auto"}. Can be set to \code{tune()} for grid search.
#' @param spectral_rank Non-negative integer rank for Spectral Local
#'   Dependence.  Zero means standard (no SLD). Can be set to
#'   \code{tune()} for grid search.
#'
#' @details
#'
#' ## What does it mean to predict?
#'
#' To predict the cluster assignment for a new observation, the model
#' computes posterior class probabilities under the estimated mixture
#' and returns the modal (most probable) class.
#'
#' ## Engine arguments
#'
#' Additional arguments may be passed to the underlying
#' \code{\link{fit_lca}()} function via
#' \code{set_engine("mixLCA", ...)}. Commonly used engine
#' arguments include:
#' \describe{
#'   \item{\code{continuous}}{Character vector of continuous indicator
#'     column names.}
#'   \item{\code{categorical}}{Character vector of categorical indicator
#'     column names.}
#'   \item{\code{concomitant}}{Character vector or one-sided formula
#'     for concomitant predictors.}
#'   \item{\code{penalty}}{Non-negative numeric or \code{"auto"}.}
#'   \item{\code{spectral_rank}}{Integer for Spectral Local
#'     Dependence rank.}
#'   \item{\code{control}}{A list from \code{\link{lca_control}()}.}
#' }
#'
#' When neither \code{continuous} nor \code{categorical} is
#' specified, column types are inferred automatically: numeric
#' columns become continuous indicators and factor or character
#' columns become categorical indicators.
#'
#' @return An \code{lca_clust} cluster specification.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("tidyclust", quietly = TRUE)) {
#'   data(voter_perceptions)
#'   set.seed(110)
#'
#'   # --- Naive LCA (local independence) ---
#'   spec_naive <- lca_clust(num_clusters = 3) |>
#'     tidyclust::set_engine("mixLCA",
#'       categorical = names(voter_perceptions),
#'       control = lca_control(n_starts = 3))
#'
#'   fit_naive <- tidyclust::fit(spec_naive, ~ ., data = voter_perceptions)
#'   predict(fit_naive, new_data = voter_perceptions)
#'
#'   # --- SLD clustering (spectral local dependence, rank 2) ---
#'   spec_sld <- lca_clust(num_clusters = 3, spectral_rank = 2L) |>
#'     tidyclust::set_engine("mixLCA",
#'       categorical = names(voter_perceptions),
#'       control = lca_control(n_starts = 3))
#'
#'   fit_sld <- tidyclust::fit(spec_sld, ~ ., data = voter_perceptions)
#'   predict(fit_sld, new_data = voter_perceptions)
#'
#'   # Compare fit indices
#'   cat("Naive BIC:", BIC(fit_naive$fit),
#'       " SLD BIC:", BIC(fit_sld$fit), "\n")
#' }
#' }
#'
#' @export
lca_clust <- function(mode = "partition",
                      engine = "mixLCA",
                      num_clusters = NULL,
                      dependence = NULL,
                      penalty = NULL,
                      spectral_rank = NULL) {
  if (!requireNamespace("tidyclust", quietly = TRUE))
    stop("Package 'tidyclust' is required for lca_clust(). ",
         "Install it with install.packages('tidyclust').",
         call. = FALSE)

  args <- list(
    num_clusters  = rlang::enquo(num_clusters),
    dependence    = rlang::enquo(dependence),
    penalty       = rlang::enquo(penalty),
    spectral_rank = rlang::enquo(spectral_rank)
  )

  tidyclust::new_cluster_spec(
    "lca_clust",
    args     = args,
    eng_args = NULL,
    mode     = mode,
    method   = NULL,
    engine   = engine
  )
}

#' @export
print.lca_clust <- function(x, ...) {
  cat("Latent Class Clustering Specification (", x$mode, ")\n\n", sep = "")
  parsnip::model_printer(x, ...)
  if (!is.null(x$method$fit$args)) {
    cat("Model fit template:\n")
    print(parsnip::show_call(x))
  }
  invisible(x)
}

# --------------------------------------------------------------------------
# Engine registration (called from .onLoad)
# --------------------------------------------------------------------------
make_lca_clust <- function() {
  modelenv::set_new_model("lca_clust")
  modelenv::set_model_mode("lca_clust", "partition")

  # -- engine --
  modelenv::set_model_engine("lca_clust", "partition", "mixLCA")
  modelenv::set_dependency(
    model = "lca_clust", mode = "partition",
    eng   = "mixLCA",    pkg  = "mixLCA"
  )

  # -- fit --
  modelenv::set_fit(
    model = "lca_clust",
    eng   = "mixLCA",
    mode  = "partition",
    value = list(
      interface = "data.frame",
      protect   = c("x", "num_clusters", "penalty", "spectral_rank"),
      func      = c(pkg = "mixLCA", fun = ".lca_clust_fit_mixLCA"),
      defaults  = list()
    )
  )

  # -- encoding --
  modelenv::set_encoding(
    model = "lca_clust",
    eng   = "mixLCA",
    mode  = "partition",
    options = list(
      predictor_indicators = "none",
      compute_intercept    = FALSE,
      remove_intercept     = TRUE,
      allow_sparse_x       = FALSE
    )
  )

  # -- model args --
  modelenv::set_model_arg(
    model    = "lca_clust",
    eng      = "mixLCA",
    exposed  = "num_clusters",
    original = "num_clusters",
    func     = list(pkg = "dials", fun = "num_clusters"),
    has_submodel = FALSE
  )
  modelenv::set_model_arg(
    model    = "lca_clust",
    eng      = "mixLCA",
    exposed  = "dependence",
    original = "dependence",
    func     = list(pkg = "base", fun = "character"),
    has_submodel = FALSE
  )
  modelenv::set_model_arg(
    model    = "lca_clust",
    eng      = "mixLCA",
    exposed  = "penalty",
    original = "penalty",
    func     = list(pkg = "dials", fun = "penalty"),
    has_submodel = FALSE
  )
  modelenv::set_model_arg(
    model    = "lca_clust",
    eng      = "mixLCA",
    exposed  = "spectral_rank",
    original = "spectral_rank",
    func     = list(pkg = "base", fun = "integer"),
    has_submodel = FALSE
  )

  # -- predict: cluster assignments --
  modelenv::set_pred(
    model = "lca_clust",
    eng   = "mixLCA",
    mode  = "partition",
    type  = "cluster",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(pkg = "mixLCA", fun = ".lca_clust_predict_mixLCA"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data)
      )
    )
  )

  # -- predict: class probabilities --
  modelenv::set_pred(
    model = "lca_clust",
    eng   = "mixLCA",
    mode  = "partition",
    type  = "prob",
    value = list(
      pre  = NULL,
      post = NULL,
      func = c(pkg = "mixLCA", fun = ".lca_clust_predict_prob_mixLCA"),
      args = list(
        object   = rlang::expr(object$fit),
        new_data = rlang::expr(new_data)
      )
    )
  )
}

# --------------------------------------------------------------------------
# .onLoad: register engine only when tidyclust is available
# --------------------------------------------------------------------------
.onLoad <- function(libname, pkgname) {
  tryCatch({
    ns_ok <- requireNamespace("modelenv", quietly = TRUE) &&
             requireNamespace("tidyclust", quietly = TRUE)
    if (ns_ok) make_lca_clust()
  }, error = function(e) invisible(NULL))
}
