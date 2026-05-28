# ==============================================================================
# File: R/01_density_continuous.R
# mixLCA - multivariate normal density with analytical missing-data
# marginalisation; M-step updates with optional L1-penalised covariance.
# ==============================================================================

#' Evaluate Continuous Log-Density
#'
#' Computes the log-density of each observation under a class-specific
#' multivariate normal distribution. When entries are missing the density
#' is evaluated over the marginal distribution of the observed subset,
#' which constitutes exact analytical marginalisation rather than
#' listwise deletion or ad hoc imputation.
#'
#' @param Y Numeric matrix (N x d), may contain \code{NA}.
#' @param mu Numeric vector of length d (class-specific mean).
#' @param Sigma Numeric d x d covariance matrix.
#' @return Numeric vector of length N (log-densities).
#' @keywords internal
eval_continuous_density <- function(Y, mu, Sigma) {
  as.vector(eval_continuous_density_cpp(Y, as.numeric(mu), as.matrix(Sigma)))
}

#' M-Step Update for Continuous Parameters
#'
#' Updates class-specific mean and covariance using posterior weights.
#' Missing values are imputed through their conditional expectations
#' under the current parameters (the standard EM treatment of
#' incomplete data). Optionally applies L1 soft-thresholding to
#' off-diagonal covariance elements, producing data-driven sparsity.
#'
#' @param Y Numeric matrix (N x d), may contain \code{NA}.
#' @param weights Numeric vector of length N (posterior weights for this
#'   class).
#' @param mu_old Current mean vector (length d).
#' @param Sigma_old Current d x d covariance matrix.
#' @param dependence Character: \code{"none"}, \code{"full"}, or
#'   \code{"penalized"}.
#' @param penalty Numeric L1 penalty (used only when
#'   \code{dependence = "penalized"}).
#' @return List with elements \code{mean} (updated mean) and
#'   \code{covariance} (updated covariance).
#' @keywords internal
update_continuous <- function(Y, weights, mu_old, Sigma_old,
                              dependence = "full", penalty = 0) {
  N <- nrow(Y)
  d <- ncol(Y)
  W <- sum(weights)

  # Accumulate sufficient statistics E[y] and E[yy']
  E_y  <- matrix(0, nrow = N, ncol = d)
  E_yy <- matrix(0, nrow = d, ncol = d)

  for (i in seq_len(N)) {
    yi  <- Y[i, ]
    mis <- which(is.na(yi))
    obs <- which(!is.na(yi))

    if (length(mis) == 0L) {
      # Complete case
      E_y[i, ] <- yi
      E_yy <- E_yy + weights[i] * tcrossprod(yi)

    } else if (length(obs) == 0L) {
      # Entirely missing: impute with current mean
      E_y[i, ] <- mu_old
      E_yy <- E_yy + weights[i] * (Sigma_old + tcrossprod(mu_old))

    } else {
      # Partial: conditional distribution of missing given observed
      S_oo <- Sigma_old[obs, obs, drop = FALSE]
      S_mo <- Sigma_old[mis, obs, drop = FALSE]
      S_mm <- Sigma_old[mis, mis, drop = FALSE]

      ridge_oo  <- max(diag(S_oo)) * 1e-6 + 1e-8
      S_oo_inv  <- solve(S_oo + diag(ridge_oo, length(obs)))
      delta     <- yi[obs] - mu_old[obs]

      cond_mean <- mu_old[mis] + S_mo %*% S_oo_inv %*% delta
      cond_cov  <- S_mm - S_mo %*% S_oo_inv %*% t(S_mo)

      completed       <- numeric(d)
      completed[obs]  <- yi[obs]
      completed[mis]  <- cond_mean
      E_y[i, ]        <- completed

      # The conditional covariance contributes to E[yy']
      C_block          <- matrix(0, d, d)
      C_block[mis, mis] <- cond_cov
      E_yy <- E_yy + weights[i] * (tcrossprod(completed) + C_block)
    }
  }

  # Updated mean and raw covariance
  new_mu    <- colSums(weights * E_y) / W
  new_Sigma <- (E_yy / W) - tcrossprod(new_mu)

  active_offdiag <- NULL

  if (dependence == "none") {
    new_Sigma <- diag(diag(new_Sigma), nrow = d)
    var_floor <- max(diag(new_Sigma)) * 1e-6 + 1e-8
    diag(new_Sigma) <- pmax(diag(new_Sigma), var_floor)

  } else if (dependence == "penalized" && penalty > 0) {
    if (!requireNamespace("glassoFast", quietly = TRUE))
      stop("Package 'glassoFast' is required for dependence = 'penalized'. ",
           "Install via install.packages('glassoFast').")
    new_Sigma <- (new_Sigma + t(new_Sigma)) / 2
    var_floor <- max(diag(new_Sigma)) * 1e-6 + 1e-8
    diag(new_Sigma) <- pmax(diag(new_Sigma), var_floor)
    gl <- glassoFast::glassoFast(new_Sigma, rho = penalty)
    new_Sigma <- gl$w
    active_offdiag <- sum(abs(gl$wi[lower.tri(gl$wi, diag = FALSE)]) > 1e-6)
    new_Sigma <- (new_Sigma + t(new_Sigma)) / 2

  } else {
    var_floor <- max(diag(new_Sigma)) * 1e-6 + 1e-8
    diag(new_Sigma) <- pmax(diag(new_Sigma), var_floor)
    new_Sigma <- (new_Sigma + t(new_Sigma)) / 2
  }

  list(mean = new_mu, covariance = new_Sigma, active_offdiag = active_offdiag)
}
