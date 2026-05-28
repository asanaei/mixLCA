# ==============================================================================
# File: R/00_utils.R
# mixLCA - core numerical utilities and input validation
# ==============================================================================

#' Numerically Stable Log-Sum-Exp
#'
#' Computes log(sum(exp(x))) without overflow or underflow.
#'
#' @param x Numeric vector.
#' @return Scalar.
#' @keywords internal
log_sum_exp <- function(x) {
  m <- max(x, na.rm = TRUE)
  if (is.infinite(m)) return(-Inf)
  m + log(sum(exp(x - m), na.rm = TRUE))
}

#' Softmax Transformation
#'
#' Row-wise softmax of a matrix. Each row sums to 1.
#'
#' @param X Numeric matrix (N x K).
#' @return Matrix of probabilities (N x K).
#' @keywords internal
softmax <- function(X) {
  row_max <- apply(X, 1, max, na.rm = TRUE)
  shifted <- exp(sweep(X, 1, row_max, "-"))
  sweep(shifted, 1, rowSums(shifted, na.rm = TRUE), "/")
}

#' Soft-Threshold Operator
#'
#' Applies element-wise soft thresholding: sign(x) * max(|x| - lambda, 0).
#' Used for L1-penalised covariance estimation in mixLCA.
#'
#' @param x Numeric value or matrix.
#' @param lambda Non-negative penalty.
#' @return Thresholded value(s).
#' @keywords internal
soft_threshold <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}

#' Count Free Parameters in a mixLCA Object
#'
#' Tallies the number of estimable parameters.
#'
#' @param model A \code{mixLCA} object.
#' @return Integer.
#' @keywords internal
count_params <- function(model) {
  K  <- model$n_classes
  np <- 0L

  # Concomitant coefficients (P predictors + intercept) x (K-1) classes
  if (!is.null(model$concomitant_coefs)) {
    np <- np + length(model$concomitant_coefs)
  } else {
    np <- np + (K - 1L)
  }

  # Continuous: class-specific means + unique covariance elements
  if (!is.null(model$specs$continuous)) {
    d <- length(model$specs$continuous)
    if (model$specs$dependence == "penalized" && !is.null(model$continuous_params)) {
      np <- np + K * d  # Class-specific means
      for (k in seq_len(K)) {
        # Use the pre-PSD active set if available; otherwise fall back to matrix inspection
        ao <- model$continuous_params$active_offdiag
        if (!is.null(ao) && !is.null(ao[[k]])) {
          np <- np + ao[[k]] + d  # cached active off-diag + diagonal variances
        } else {
          S <- model$continuous_params$covariances[[k]]
          np <- np + sum(abs(S[lower.tri(S, diag = FALSE)]) > 1e-4) + d
        }
      }
    } else {
      cov_per_class <- if (model$specs$dependence == "none") d else d * (d + 1L) / 2L
      np <- np + K * (d + cov_per_class)
    }
  }

  # Categorical: (n_categories - 1) per variable per class
  if (!is.null(model$categorical_params)) {
    for (k in seq_len(K)) {
      for (vname in names(model$categorical_params[[k]])) {
        np <- np + length(model$categorical_params[[k]][[vname]]) - 1L
      }
    }
  }

  # Direct effects: additional (C_parent - 1) * (C_child - 1) per class per pair
  if (!is.null(model$cat_dep_params) && !is.null(model$specs$cat_direct_effects)) {
    for (d in seq_along(model$specs$cat_direct_effects)) {
      # cond_mat from class 1 gives us category counts
      cond_mat <- model$cat_dep_params[[1]][[d]]
      C_parent <- nrow(cond_mat)
      C_child  <- ncol(cond_mat)
      # Each parent category gives (C_child - 1) free params;
      # the marginal (C_child - 1) is already counted above.
      # Net addition = (C_parent - 1) * (C_child - 1) per class.
      np <- np + K * (C_parent - 1L) * (C_child - 1L)
    }
  }

  # Spectral Local Dependence Grassmannian parameter count
  if (isTRUE(any(model$specs$spectral_rank > 0L)) && !is.null(model$categorical_params)) {
    d_vec <- as.integer(model$specs$spectral_rank)
    if (length(d_vec) == 1L) d_vec <- rep(d_vec, model$n_classes)

    J <- length(model$specs$categorical)
    C <- sum(sapply(model$categorical_params[[1]], length))

    if (isTRUE(model$specs$spectral_pool)) {
      d <- max(d_vec)
      d_eff <- min(d, max(C - J, 0L))
      basis_df <- max(d_eff * (C - J - d_eff), 0L)
      if (basis_df > 0L) np <- np + basis_df
    } else {
      for (k in seq_len(model$n_classes)) {
        d <- d_vec[k]
        d_eff <- min(d, max(C - J, 0L))
        basis_df <- max(d_eff * (C - J - d_eff), 0L)
        if (basis_df > 0L) np <- np + basis_df
      }
    }
  }

  as.integer(np)
}

#' Validate Inputs for lca
#'
#' Checks data, variable names, and parameter sanity before estimation.
#'
#' @param data Data frame.
#' @param continuous Character vector or NULL.
#' @param categorical Character vector or NULL.
#' @param concomitant Character vector or NULL.
#' @param n_classes Positive integer.
#' @param dependence One of \code{"none"}, \code{"full"}, \code{"penalized"}.
#' @param penalty Non-negative numeric.
#' @return Invisible TRUE on success; stops with an informative message
#'   otherwise.
#' @keywords internal
validate_inputs <- function(data, continuous, categorical, concomitant,
                            n_classes, dependence, penalty) {

  if (!is.data.frame(data))
    stop("`data` must be a data frame.")
  if (nrow(data) < 10L)
    stop("Insufficient observations (minimum 10).")
  if (is.null(continuous) && is.null(categorical))
    stop("Specify at least one of `continuous` or `categorical` indicators.")

  concom_vars <- if (inherits(concomitant, "formula"))
    all.vars(concomitant) else concomitant

  all_vars     <- c(continuous, categorical, concom_vars)
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0L)
    stop("Variables absent from data: ", paste(missing_vars, collapse = ", "))

  if (!is.null(concomitant)) {
    if (inherits(concomitant, "formula")) {
      mf <- tryCatch(
        stats::model.frame(concomitant, data = data,
                           na.action = stats::na.pass),
        error = function(e) NULL)
      has_na <- !is.null(mf) && any(!stats::complete.cases(mf))
    } else {
      has_na <- any(!stats::complete.cases(data[, concomitant, drop = FALSE]))
    }
    if (has_na)
      stop("Missing values detected in concomitant predictors. ",
           "Impute or filter these rows before calling fit_lca() so that ",
           "row alignment between data and posteriors is preserved.",
           call. = FALSE)
  }

  if (!is.numeric(n_classes) || length(n_classes) != 1L || n_classes < 1L)
    stop("`n_classes` must be an integer >= 1.")

  valid_dep <- c("none", "full", "penalized")
  if (!(dependence %in% valid_dep))
    stop("`dependence` must be one of: ", paste(valid_dep, collapse = ", "))

  if (dependence == "penalized" && (!is.numeric(penalty) || penalty < 0))
    stop("`penalty` must be a non-negative number when dependence = 'penalized'.")

  if (!is.null(continuous)) {
    for (v in continuous)
      if (!is.numeric(data[[v]]))
        stop("Continuous indicator '", v, "' is not numeric.")
  }
  if (!is.null(categorical)) {
    for (v in categorical)
      if (length(unique(stats::na.omit(data[[v]]))) < 2L)
        stop("Categorical indicator '", v, "' has fewer than 2 levels.")
  }

  invisible(TRUE)
}
if(getRversion() >= "2.15.1") utils::globalVariables(c("magnitude", "name", "iteration", "log_lik", "value", "weight", "indicator", "lower", "upper", "max_posterior", "Category", "Probability", "Class", "Variable", "loading", "column_label", "item", "index", "eigenvalue", "class"))
