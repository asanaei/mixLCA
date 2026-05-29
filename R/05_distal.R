# ==============================================================================
# File: R/05_distal.R
# mixLCA - distal outcome estimation via the BCH method.
#
# The classification error matrix is computed from frozen posteriors and
# then inverted to yield corrective observation weights. No gradient from
# distal outcomes flows back into the measurement model (the "cut").
# ==============================================================================

#' Compute BCH Classification Error Matrix
#'
#' Constructs the K x K matrix E where E[w, k] = P(W = w | C = k),
#' i.e. the probability that an observation truly belonging to class k
#' is modally assigned to class w.
#'
#' @param posteriors N x K posterior probability matrix.
#' @return K x K classification error matrix.
#' @keywords internal
classification_error_matrix <- function(posteriors) {
  K     <- ncol(posteriors)
  modal <- apply(posteriors, 1, which.max)
  N_k   <- colSums(posteriors)
  E     <- matrix(0, K, K)

  for (w in seq_len(K)) {
    idx <- which(modal == w)
    if (length(idx) > 0L) {
      for (k in seq_len(K)) {
        E[w, k] <- sum(posteriors[idx, k]) / N_k[k]
      }
    }
  }
  # Safety: ensure rows with zero modal assignment get identity
  for (w in seq_len(K)) {
    if (sum(E[w, ]) < 1e-15) E[w, w] <- 1
  }
  E
}

#' Compute BCH Observation Weights
#'
#' Inverts the classification error matrix and produces per-observation,
#' per-class weights that correct for classification uncertainty.
#'
#' @param posteriors N x K posterior probability matrix.
#' @return List with elements \code{weights} (N x K), \code{W_inv}
#'   (K x K inverse error matrix), and \code{error_matrix} (K x K).
#' @keywords internal
bch_weights <- function(posteriors) {
  K     <- ncol(posteriors)
  N     <- nrow(posteriors)
  E     <- classification_error_matrix(posteriors)
  W_inv <- solve(E + diag(1e-6, K))
  modal <- apply(posteriors, 1, which.max)

  wt <- matrix(0, N, K)
  for (i in seq_len(N))
    wt[i, ] <- W_inv[, modal[i]]

  list(weights = wt, W_inv = W_inv, error_matrix = E)
}

#' Distal Outcome Estimation via BCH Weighting
#'
#' Estimates class-specific regression models for one or more distal
#' outcomes, using BCH inverse-classification-error weights. The
#' measurement model posteriors are fixed before this step (the "cut"):
#' no gradient from distal outcomes reaches the class definitions.
#'
#' @param model Fitted \code{mixLCA} object.
#' @param data Data frame containing both the original variables and the
#'   distal outcome.
#' @param formula A formula for the distal model, e.g.\cr
#'   \code{outcome ~ predictor1 + predictor2}.\cr
#'   A right-hand side of \code{~ 1} estimates unconditional class
#'   means.
#' @param family Character: \code{"gaussian"}, \code{"binomial"}, or
#'   \code{"poisson"}.
#' @return An object of class \code{mixDistal} containing class-specific
#'   model summaries.
#'
#' @references
#' Bolck, A., Croon, M., & Hagenaars, J. (2004). Estimating latent
#' structure models with categorical variables: One-step versus
#' three-step estimators. \emph{Political Analysis}, 12(1), 3-27.
#'
#' Vermunt, J. K. (2010). Latent class modeling with covariates: Two
#' improved three-step approaches. \emph{Political Analysis}, 18(4),
#' 450-469.
#'
#' Bakk, Z., Tekle, F. B., & Vermunt, J. K. (2013). Estimating the
#' association between latent class membership and external variables
#' using bias-adjusted three-step approaches. \emph{Sociological
#' Methodology}, 43(1), 272-311.
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
#'
#' # Binary distal outcome under BCH weighting
#' d_bin <- distal(fit, health_screening,
#'                 outcome ~ age, family = "binomial")
#' print(d_bin)
#'
#' # Continuous distal outcome (uses marker_4 as a stand-in)
#' d_gauss <- distal(fit, health_screening,
#'                   marker_4 ~ age, family = "gaussian")
#' print(d_gauss)
#' }
#' @export
distal <- function(model, data, formula, family = "gaussian") {
  if (!inherits(model, "mixLCA"))
    stop("`model` must be a mixLCA object.")

  posteriors <- model$posteriors
  K  <- model$n_classes
  N  <- nrow(posteriors)

  bch <- bch_weights(posteriors)
  wt  <- bch$weights

  # Parse formula
  mf         <- stats::model.frame(formula, data = data,
                                   na.action = stats::na.pass)
  response   <- stats::model.response(mf)
  predictors <- stats::model.matrix(formula, data = mf)

  # Identify complete cases
  complete <- stats::complete.cases(cbind(response, predictors))

  # For binomial, coerce a factor / logical response to {0,1}: take
  # the second level as the "success" (matches stats::glm convention).
  if (family == "binomial" && is.factor(response)) {
    if (nlevels(response) != 2L)
      stop("binomial distal response must be a two-level factor.")
    response <- as.integer(response) - 1L
  } else if (family == "binomial" && is.logical(response)) {
    response <- as.integer(response)
  }

  results <- list()
  for (k in seq_len(K)) {
    w_k <- wt[complete, k]
    Xc  <- predictors[complete, , drop = FALSE]
    yc  <- response[complete]

    if (family == "gaussian") {
      # Manual WLS: BCH weights may be negative, so lm.wfit is unsuitable
      XtWX     <- crossprod(Xc * w_k, Xc)
      XtWy     <- crossprod(Xc * w_k, yc)
      XtWX_inv <- MASS::ginv(XtWX)

      coef_est  <- as.vector(XtWX_inv %*% XtWy)
      names(coef_est) <- colnames(Xc)
      residuals <- yc - Xc %*% coef_est
      W_sum     <- sum(w_k)
      sigma2    <- sum(w_k * residuals^2) / W_sum

      # Sandwich standard errors
      scores    <- Xc * as.vector(w_k * residuals)
      meat      <- crossprod(scores)
      vcov_sand <- XtWX_inv %*% meat %*% XtWX_inv
      var_diag  <- diag(vcov_sand)
      if (any(var_diag < 0)) {
        warning(sprintf(
          "Class %d (gaussian): negative variance(s) in BCH sandwich estimator; standard errors set to NA. The distal model may be unstable for this class.",
          k), call. = FALSE)
        se <- rep(NA_real_, length(var_diag))
      } else {
        se <- sqrt(var_diag)
      }
      z_val <- ifelse(is.na(se), NA_real_, coef_est / se)
      p_val <- ifelse(is.na(se), NA_real_, 2 * stats::pnorm(-abs(z_val)))

      coef_table <- cbind(Estimate = coef_est, SE = se,
                          z = z_val, p = p_val)

      results[[k]] <- list(
        class      = k,
        family     = family,
        coef_table = coef_table,
        sigma      = sqrt(abs(sigma2)),
        n_eff      = W_sum
      )

    } else {
      # BCH weights can be negative. Standard optim/glm diverges with negative
      # weights because the log-likelihood becomes unbounded. We solve the
      # weighted score equations directly via IRLS.
      p_dim <- ncol(Xc)
      coef_est <- rep(0, p_dim)
      irls_converged <- FALSE

      for (irls_iter in seq_len(50L)) {
        eta <- as.vector(Xc %*% coef_est)

        if (family == "binomial") {
          mu_val  <- 1 / (1 + exp(-eta))
          var_val <- mu_val * (1 - mu_val)
        } else if (family == "poisson") {
          mu_val  <- exp(pmin(eta, 20))
          var_val <- mu_val
        } else {
          stop("Unsupported family: ", family)
        }

        var_val <- pmax(var_val, 1e-10)
        resid   <- yc - mu_val
        score   <- crossprod(Xc, w_k * resid)

        W_diag <- w_k * var_val
        XtWX   <- crossprod(Xc * W_diag, Xc)

        # Project XtWX to positive definite (negative BCH weights can flip eigenvalues)
        eig_X <- eigen(XtWX, symmetric = TRUE)
        eig_ridge <- max(abs(eig_X$values)) * 1e-6 + 1e-8
        eig_X$values <- pmax(abs(eig_X$values), eig_ridge)
        XtWX_pd <- eig_X$vectors %*% diag(eig_X$values, nrow = length(eig_X$values)) %*% t(eig_X$vectors)
        delta <- as.vector(solve(XtWX_pd, score))

        # Step-halving: cut step until the next linear predictor stays bounded
        step_size <- 1.0
        for (half in seq_len(5L)) {
          if (max(abs(Xc %*% (coef_est + step_size * delta))) < 30) break
          step_size <- step_size * 0.5
        }
        coef_est <- coef_est + step_size * delta
        if (max(abs(step_size * delta)) < 1e-8) { irls_converged <- TRUE; break }
      }

      if (!irls_converged) {
        results[[k]] <- list(
          class = k, family = family, coef_table = NULL,
          note  = "IRLS did not converge for BCH-weighted GLM."
        )
        next
      }

      names(coef_est) <- colnames(Xc)
      eta <- as.vector(Xc %*% coef_est)
      if (family == "binomial") {
        mu_val  <- 1 / (1 + exp(-eta))
        var_val <- mu_val * (1 - mu_val)
      } else {
        mu_val  <- exp(pmin(eta, 20))
        var_val <- mu_val
      }
      var_val <- pmax(var_val, 1e-10)

      # Sandwich standard errors
      W_diag   <- w_k * var_val
      XtWX_raw <- crossprod(Xc * W_diag, Xc)
      eig_X    <- eigen(XtWX_raw, symmetric = TRUE)
      eig_ridge <- max(abs(eig_X$values)) * 1e-6 + 1e-8
      eig_X$values <- pmax(abs(eig_X$values), eig_ridge)
      XtWX_pd  <- eig_X$vectors %*% diag(eig_X$values, nrow = length(eig_X$values)) %*% t(eig_X$vectors)
      XtWX_inv <- solve(XtWX_pd)
      resid_dev <- yc - mu_val
      scores    <- Xc * as.vector(w_k * resid_dev)
      meat      <- crossprod(scores)
      vcov_sand <- XtWX_inv %*% meat %*% XtWX_inv
      var_diag  <- diag(vcov_sand)
      if (any(var_diag < 0)) {
        warning(sprintf(
          "Class %d (%s): negative variance(s) in BCH IRLS sandwich estimator; standard errors set to NA. The distal model may be unstable for this class.",
          k, family), call. = FALSE)
        se <- rep(NA_real_, p_dim)
      } else {
        se <- sqrt(var_diag)
      }

      z_val      <- ifelse(is.na(se), NA_real_, coef_est / se)
      p_val      <- ifelse(is.na(se), NA_real_, 2 * stats::pnorm(-abs(z_val)))
      coef_table <- cbind(Estimate = coef_est, SE = se,
                          z = z_val, p = p_val)

      results[[k]] <- list(
        class      = k,
        family     = family,
        coef_table = coef_table,
        n_eff      = sum(w_k)
      )
    }
  }

  out <- list(
    class_models = results,
    n_classes    = K,
    formula      = formula,
    family       = family,
    error_matrix = bch$error_matrix,
    W_inv        = bch$W_inv
  )
  class(out) <- "mixDistal"
  out
}

#' Print a mixDistal Object
#'
#' Prints the per-class regression coefficients, sandwich standard
#' errors, and classification error matrix produced by
#' \code{\link{distal}}.
#'
#' @param x A \code{mixDistal} object.
#' @param ... Unused.
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' \donttest{
#' data(health_screening)
#' fit <- fit_lca(health_screening,
#'                continuous  = c("marker_1","marker_2","marker_3","marker_4"),
#'                concomitant = ~ age,
#'                n_classes   = 2,
#'                control     = lca_control(n_starts = 2),
#'                verbose     = FALSE)
#' d <- distal(fit, health_screening, outcome ~ age, family = "binomial")
#' print(d)
#' }
#' @export
print.mixDistal <- function(x, ...) {
  cat("\nDistal Outcome Estimation (BCH Method) - mixLCA\n")
  cat("================================================\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Family :", x$family, "\n")
  cat("Classes:", x$n_classes, "\n\n")

  for (k in seq_len(x$n_classes)) {
    m <- x$class_models[[k]]
    cat("--- Class", k, "---\n")
    if (is.null(m$coef_table)) {
      cat("  ", m$note, "\n\n")
    } else {
      stats::printCoefmat(m$coef_table, digits = 4,
                          signif.stars = TRUE,
                          P.values = TRUE, has.Pvalue = TRUE)
      if (!is.null(m$sigma))
        cat("  Residual SD:", round(m$sigma, 4), "\n")
      cat("  Effective N:", round(m$n_eff, 1), "\n\n")
    }
  }

  cat("Classification Error Matrix:\n")
  print(round(x$error_matrix, 4))
  cat("\n")
  invisible(x)
}
