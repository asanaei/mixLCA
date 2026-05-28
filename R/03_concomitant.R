# ==============================================================================
# File: R/03_concomitant.R
# mixLCA - multinomial logistic regression for class membership probabilities.
# Reference class is always class 1 (gamma_1 = 0 by construction).
# ==============================================================================

#' Negative Log-Likelihood for Concomitant Model
#'
#' @param par Numeric vector of length P*(K-1), stored column-major.
#' @param X Numeric matrix (N x P), includes intercept column.
#' @param posteriors Numeric matrix (N x K) of current posterior
#'   probabilities.
#' @return Scalar (negative expected log-likelihood contribution).
#' @keywords internal
concomitant_nll <- function(par, X, posteriors) {
  K    <- ncol(posteriors)
  P    <- ncol(X)
  B    <- matrix(par, nrow = P, ncol = K - 1L)
  eta  <- cbind(0, X %*% B)
  pi_mat <- softmax(eta)
  -sum(posteriors * log(pi_mat + 1e-300)) + 1e-4 * sum(B^2)
}

#' Analytical Gradient for Concomitant Model
#'
#' @param par Numeric vector of length P*(K-1).
#' @param X Numeric matrix (N x P).
#' @param posteriors Numeric matrix (N x K).
#' @return Numeric vector (gradient of the negative log-likelihood).
#' @keywords internal
concomitant_grad <- function(par, X, posteriors) {
  K    <- ncol(posteriors)
  P    <- ncol(X)
  B    <- matrix(par, nrow = P, ncol = K - 1L)
  eta  <- cbind(0, X %*% B)
  pi_mat <- softmax(eta)

  G <- matrix(0, nrow = P, ncol = K - 1L)
  for (k in seq_len(K - 1L)) {
    residual <- posteriors[, k + 1L] - pi_mat[, k + 1L]
    G[, k]   <- -crossprod(X, residual) + 2e-4 * B[, k]
  }
  as.vector(G)
}

#' M-Step: Update Concomitant Coefficients
#'
#' Wraps BFGS optimisation of the multinomial logistic model.
#'
#' @param X Numeric matrix (N x P), includes intercept.
#' @param posteriors Numeric matrix (N x K).
#' @param current_coefs Current P x (K-1) coefficient matrix.
#' @return Updated P x (K-1) coefficient matrix.
#' @keywords internal
update_concomitant <- function(X, posteriors, current_coefs) {
  opt <- stats::optim(
    par      = as.vector(current_coefs),
    fn       = concomitant_nll,
    gr       = concomitant_grad,
    X        = X,
    posteriors = posteriors,
    method   = "BFGS",
    control  = list(maxit = 20L)
  )
  out <- matrix(opt$par, nrow = ncol(X), ncol = ncol(posteriors) - 1L)
  rownames(out) <- colnames(X)
  out
}

#' Compute Prior Probabilities from the Concomitant Model
#'
#' @param X Numeric matrix (N x P), includes intercept.
#' @param coefs P x (K-1) coefficient matrix.
#' @return Numeric matrix (N x K) of prior class probabilities.
#' @keywords internal
compute_priors <- function(X, coefs) {
  eta <- cbind(0, X %*% coefs)
  softmax(eta)
}
