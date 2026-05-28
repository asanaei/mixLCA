# ==============================================================================
# File: R/02_density_categorical.R
# mixLCA - categorical (product-multinomial) density and M-step updates.
# Missing categories contribute zero to the log-density (marginalisation).
# ==============================================================================

#' Evaluate Categorical Log-Density
#'
#' Computes the log-probability of each observation's categorical
#' indicators under class-specific multinomial parameters. Missing
#' values are marginalised out (contribute zero to the sum).
#'
#' @param df Data frame of categorical indicators (N x J).
#' @param probs Named list: each element is a named probability vector
#'   for one variable.
#' @return Numeric vector of length N (log-densities).
#' @keywords internal
eval_categorical_density <- function(df, probs) {
  N        <- nrow(df)
  J        <- ncol(df)
  varnames <- colnames(df)
  ld       <- numeric(N)

  for (j in seq_len(J)) {
    vname  <- varnames[j]
    values <- as.character(df[[j]])
    p_vec  <- probs[[vname]]

    assigned            <- p_vec[values]
    assigned[is.na(assigned)] <- NA_real_
    assigned[!is.na(assigned) & assigned < 1e-15] <- 1e-15

    contrib <- log(assigned)
    contrib[is.na(values) | is.na(contrib)] <- 0
    ld <- ld + contrib
  }

  ld
}

#' M-Step Update for Categorical Parameters
#'
#' Re-estimates class-conditional category probabilities using posterior
#' weights. Missing values are excluded from both numerator and
#' denominator.
#'
#' @param df Data frame of categorical indicators (N x J).
#' @param weights Numeric vector of length N (posterior weights for this
#'   class).
#' @return Named list of named probability vectors (one per variable).
#' @keywords internal
update_categorical <- function(df, weights) {
  J        <- ncol(df)
  varnames <- colnames(df)
  out      <- list()

  for (j in seq_len(J)) {
    vname   <- varnames[j]
    values  <- as.character(df[[j]])
    cats    <- sort(unique(stats::na.omit(values)))

    valid   <- !is.na(values)
    W_valid <- sum(weights[valid])

    p_vec       <- numeric(length(cats))
    names(p_vec) <- cats

    if (W_valid < 1e-15) {
      warning(sprintf("Class has near-zero valid weight for categorical variable '%s'; defaulting to uniform probabilities.", vname), call. = FALSE)
      p_vec[] <- 1.0 / length(cats)
    } else {
      for (cc in cats) {
        hit      <- valid & (values == cc)
        p_vec[cc] <- sum(weights[hit]) / W_valid
      }
    }

    # Numerical floor and renormalise
    p_vec        <- pmax(p_vec, 1e-10)
    p_vec        <- p_vec / sum(p_vec)
    out[[vname]] <- p_vec
  }

  out
}
