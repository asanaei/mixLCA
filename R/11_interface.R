# ==============================================================================
# File: R/11_interface.R
# mixLCA - user-facing entry point.
# ==============================================================================

#' Control Parameters for mixLCA
#'
#' Bundles optimiser settings into a single list for use with
#' \code{\link{fit_lca}}. Defaults reproduce the behaviour of prior
#' beta releases.
#'
#' @param max_iter Maximum EM iterations per start.
#' @param tol Convergence tolerance on absolute log-likelihood change.
#' @param n_starts Integer: number of random starting points. For
#'   publication-quality results consider at least 10.
#' @param seed Base random seed; start \emph{s} uses \code{seed + s}.
#'   The global \code{.Random.seed} is never modified.
#' @param kmeans_nstart Integer: random initializations for the internal
#'   \code{kmeans} used to seed continuous-indicator starting values.
#'   Irrelevant for categorical-only models or when \code{init_model}
#'   is supplied.
#' @return A list of control values.
#' @export
lca_control <- function(max_iter = 500L, tol = 1e-6,
                        n_starts = 1L, seed = 110L,
                        kmeans_nstart = 1L) {
  list(max_iter = as.integer(max_iter),
       tol = as.numeric(tol),
       n_starts = as.integer(n_starts),
       seed = as.integer(seed),
       kmeans_nstart = as.integer(kmeans_nstart))
}

#' Fit a Latent Class Model with mixLCA
#'
#' Estimates a finite mixture model with a partitioned architecture.
#' Antecedent concomitant predictors enter via multinomial logistic
#' regression; contemporaneous manifest indicators (continuous,
#' categorical, or mixed) define the measurement model; subsequent
#' distal outcomes are estimated separately via \code{\link{distal}}
#' under BCH inverse-classification-error weighting.
#'
#' The measurement likelihood is maximised by EM. Continuous indicators
#' follow class-specific multivariate normal distributions; missing
#' entries are marginalised analytically. Categorical indicators follow
#' a product-multinomial distribution; missing categories are omitted
#' from the likelihood contribution. Non-diagonal covariance matrices
#' may be freely estimated (\code{dependence = "full"}) or estimated
#' with sparsity through the graphical lasso
#' (\code{dependence = "penalized"}).
#'
#' Parallel execution of multiple starts is delegated to the user's
#' active \code{future::plan()}. Set the plan in your session
#' (\code{future::plan(future::multisession, workers = ...)}) before
#' calling \code{fit_lca}; otherwise starts run sequentially.
#'
#' @param data Data frame of observations.
#' @param continuous Character vector of continuous indicator variable
#'   names, or NULL.
#' @param categorical Character vector of categorical indicator variable
#'   names, or NULL.
#' @param concomitant Character vector of concomitant predictor names,
#'   or a one-sided formula (e.g.\ \code{~ age * income + poly(bmi, 2)}),
#'   or NULL. Missing values in concomitant predictors are not
#'   permitted; impute or filter rows before calling.
#' @param n_classes Integer >= 2: number of latent classes.
#' @param dependence Character controlling the covariance structure of
#'   continuous indicators within each class:
#'   \describe{
#'     \item{\code{"none"}}{Diagonal covariance (local independence).}
#'     \item{\code{"full"}}{Unrestricted covariance.}
#'     \item{\code{"penalized"}}{Graphical-lasso penalised covariance,
#'       guaranteeing exact sparsity and positive definiteness.}
#'   }
#' @param penalty Penalty for \code{dependence = "penalized"}. Either a
#'   non-negative scalar or the string \code{"auto"} (default). When
#'   \code{"auto"}, a heuristic value is selected from the data; when
#'   numeric, the supplied value is respected exactly (including 0,
#'   which then yields no shrinkage).
#' @param control List of optimiser settings; see \code{\link{lca_control}}.
#' @param verbose Logical: print per-start progress?
#' @param cat_direct_effects List of two-element character vectors
#'   specifying direct effects between categorical indicators to
#'   address local dependence violations. Each element
#'   \code{c("parent", "child")} allows the child's response
#'   probabilities to depend on the parent's observed value within
#'   each class (Vermunt, 1999). Use \code{\link{bvr_categorical}}
#'   to identify candidate pairs. When \code{NULL} (default), standard
#'   local independence is assumed for categorical indicators.
#' @param spectral_rank Integer (scalar or length-K vector): target
#'   rank \emph{d} for Spectral Local Dependence (SLD) among
#'   categorical items. If any element is > 0, SLD is activated.
#' @param spectral_pool Logical: Pool the conditional Burt matrices
#'   across classes.
#' @param init_model Optional \code{mixLCA} object for warm-starting.
#'   The prior model must have the same \code{n_classes}.
#'
#' @return An opaque object of class \code{mixLCA}. Use the provided
#'   accessors and S3 methods (\code{coef()}, \code{predict()},
#'   \code{summary()}, \code{get_posteriors()}, \code{get_loadings()},
#'   \code{plot()}) rather than \code{$}-indexing internal fields,
#'   which may be restructured in future versions. Key downstream
#'   functions: \code{\link{distal}}, \code{\link{fit_indices}},
#'   \code{\link{bvr_tests}}, \code{\link{enumerate_lca}}.
#'
#' @examples
#' \dontrun{
#' withr::with_seed(110, {
#'   N  <- 600
#'   cl <- sample(1:2, N, replace = TRUE)
#'   df <- data.frame(
#'     x1  = ifelse(cl == 1, rnorm(N, 10, 2), rnorm(N, 5, 2)),
#'     x2  = ifelse(cl == 1, rnorm(N, 10, 2), rnorm(N, 5, 2)),
#'     x3  = ifelse(cl == 1, rnorm(N,  3, 1), rnorm(N, 8, 1)),
#'     age = rnorm(N, 40, 10),
#'     outcome = ifelse(cl == 1, rnorm(N, 80, 10), rnorm(N, 40, 10))
#'   )
#' })
#'
#' fit <- fit_lca(
#'   data        = df,
#'   continuous  = c("x1", "x2", "x3"),
#'   concomitant = ~ age,
#'   n_classes   = 2,
#'   dependence  = "full",
#'   control     = lca_control(n_starts = 5)
#' )
#'
#' summary(fit, data = df)
#' fit_indices(fit)
#' bvr_tests(fit, df)
#' plot(fit, type = "profiles")
#' plot(fit, type = "bvr", data = df)
#'
#' dis <- distal(fit, df, outcome ~ age, family = "gaussian")
#' print(dis)
#' }
#'
#' @export
fit_lca <- function(data,
                    continuous  = NULL,
                    categorical = NULL,
                    concomitant = NULL,
                    n_classes   = 2L,
                    dependence  = "full",
                    penalty     = "auto",
                    control     = lca_control(),
                    verbose     = TRUE,
                    cat_direct_effects = NULL,
                    spectral_rank = 0L,
                    spectral_pool = FALSE,
                    init_model  = NULL) {

  cl        <- match.call()
  n_classes <- as.integer(n_classes)

  if (!is.list(control) ||
      !all(c("max_iter", "tol", "n_starts", "seed",
             "kmeans_nstart") %in% names(control)))
    stop("`control` must be a list produced by lca_control().")

  # Auto-select penalty
  if (dependence == "penalized" && identical(penalty, "auto")) {
    N <- nrow(data)
    if (!is.null(continuous)) {
      vars    <- sapply(continuous, function(v)
        stats::var(data[[v]], na.rm = TRUE))
      penalty <- median(vars) / sqrt(N)
    } else {
      penalty <- 0.01
    }
    if (verbose)
      message("Auto-selected penalty = ", round(penalty, 5))
  } else if (identical(penalty, "auto")) {
    # auto requested but dependence != penalized; collapse to numeric 0
    penalty <- 0
  }

  if (!is.numeric(penalty) || length(penalty) != 1L || penalty < 0)
    stop("`penalty` must be a non-negative scalar or the string \"auto\".")

  validate_inputs(data, continuous, categorical, concomitant,
                  n_classes, dependence, penalty)

  user_set_dep <- "dependence" %in% names(cl) || "penalty" %in% names(cl)
  if (user_set_dep && is.null(continuous) &&
      (dependence != "none" || penalty > 0)) {
    message("Note: `dependence` and `penalty` only affect continuous ",
            "indicators. With categorical-only models they have no effect.")
  }

  if (is.numeric(spectral_rank) && any(spectral_rank > 0L)) {
    if (is.null(categorical) || length(categorical) < 2L)
      stop("Spectral local dependence requires at least two categorical indicators.")
    if (!is.null(cat_direct_effects)) {
      warning("Both 'spectral_rank' and 'cat_direct_effects' are specified. ",
              "Spectral local dependence supersedes pairwise direct effects. ",
              "Ignoring direct effects.")
      cat_direct_effects <- NULL
    }
  }

  model <- if (control$n_starts == 1L) {
    run_em(
      data        = data,
      continuous  = continuous,
      categorical = categorical,
      concomitant = concomitant,
      n_classes   = n_classes,
      dependence  = dependence,
      penalty     = penalty,
      max_iter    = control$max_iter,
      tol         = control$tol,
      seed        = control$seed,
      cat_direct_effects = cat_direct_effects,
      spectral_rank = spectral_rank,
      spectral_pool = spectral_pool,
      init_model  = init_model,
      kmeans_nstart = control$kmeans_nstart
    )
  } else {
    run_em_robust(
      data        = data,
      continuous  = continuous,
      categorical = categorical,
      concomitant = concomitant,
      n_classes   = n_classes,
      dependence  = dependence,
      penalty     = penalty,
      max_iter    = control$max_iter,
      tol         = control$tol,
      n_starts    = control$n_starts,
      base_seed   = control$seed,
      verbose     = verbose,
      cat_direct_effects = cat_direct_effects,
      spectral_rank = spectral_rank,
      spectral_pool = spectral_pool,
      init_model  = init_model,
      kmeans_nstart = control$kmeans_nstart
    )
  }

  model$call <- cl
  model
}
