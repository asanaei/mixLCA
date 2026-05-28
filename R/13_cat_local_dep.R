# ==============================================================================
# File: R/13_cat_local_dep.R
# mixLCA - Categorical local dependence via direct effects
# (Hagenaars-Vermunt-Magidson specification search)
# ==============================================================================

#' Bivariate Residuals for Categorical Indicators
#'
#' Computes the Pearson chi-squared bivariate residual for each pair of
#' categorical manifest variables. Under local independence the model-implied
#' bivariate frequencies should match the observed ones. Large BVR values
#' (rule of thumb > 4) signal local dependence addressable by adding direct
#' effects via \code{cat_direct_effects} in \code{\link{fit_lca}}.
#'
#' @param model A \code{mixLCA} object.
#' @param data Data frame used for estimation.
#' @return Data frame with columns \code{var1}, \code{var2}, \code{bvr},
#'   \code{df}, \code{p_value}, ordered by descending BVR.
#'
#' @examples
#' \donttest{
#' data(voter_perceptions)
#' fit <- fit_lca(voter_perceptions,
#'                categorical = names(voter_perceptions),
#'                n_classes   = 3,
#'                control     = lca_control(n_starts = 2, seed = 110),
#'                verbose     = FALSE)
#' head(bvr_categorical(fit, voter_perceptions), 5)
#' }
#' @export
bvr_categorical <- function(model, data) {
  if (!inherits(model, "mixLCA"))
    stop("`model` must be a mixLCA object.")

  if (isTRUE(any(model$specs$spectral_rank > 0L))) {
    warning("Model includes spectral local dependence. BVR expected frequencies reflect residuals from the unshifted marginal profiles.")
  }

  cat_vars <- model$specs$categorical
  if (is.null(cat_vars) || length(cat_vars) < 2)
    stop("At least 2 categorical indicators are required.")

  K  <- model$n_classes
  pi_k <- colMeans(model$posteriors)
  cat_params <- model$categorical_params

  pairs <- utils::combn(cat_vars, 2, simplify = FALSE)
  rows  <- vector("list", length(pairs))

  for (idx in seq_along(pairs)) {
    v1 <- pairs[[idx]][1]
    v2 <- pairs[[idx]][2]

    vals1 <- as.character(data[[v1]])
    vals2 <- as.character(data[[v2]])

    ok    <- !is.na(vals1) & !is.na(vals2)
    v1c   <- vals1[ok]
    v2c   <- vals2[ok]
    n_ok  <- sum(ok)

    cats1 <- sort(unique(v1c))
    cats2 <- sort(unique(v2c))

    O <- table(factor(v1c, levels = cats1), factor(v2c, levels = cats2))

    E <- matrix(0, length(cats1), length(cats2),
                dimnames = list(cats1, cats2))

    # Check if this pair already has a direct effect modeled
    has_dir <- FALSE
    dir_idx <- 0L
    is_v1_parent <- FALSE
    if (!is.null(model$specs$cat_direct_effects)) {
      for (d in seq_along(model$specs$cat_direct_effects)) {
        pair <- model$specs$cat_direct_effects[[d]]
        if (pair[1] == v1 && pair[2] == v2) { has_dir <- TRUE; dir_idx <- d; is_v1_parent <- TRUE; break }
        if (pair[1] == v2 && pair[2] == v1) { has_dir <- TRUE; dir_idx <- d; is_v1_parent <- FALSE; break }
      }
    }

    for (k in seq_len(K)) {
      p1 <- cat_params[[k]][[v1]]
      p2 <- cat_params[[k]][[v2]]
      for (a in cats1) {
        for (b in cats2) {
          if (has_dir && !is.null(model$cat_dep_params)) {
            cond_mat <- model$cat_dep_params[[k]][[dir_idx]]
            if (is_v1_parent) {
              pa <- if (a %in% names(p1)) p1[[a]] else 0
              pb <- if (a %in% rownames(cond_mat) && b %in% colnames(cond_mat)) cond_mat[a, b] else 0
            } else {
              pb <- if (b %in% names(p2)) p2[[b]] else 0
              pa <- if (b %in% rownames(cond_mat) && a %in% colnames(cond_mat)) cond_mat[b, a] else 0
            }
            E[a, b] <- E[a, b] + pi_k[k] * pa * pb
          } else {
            pa <- if (a %in% names(p1)) p1[[a]] else 0
            pb <- if (b %in% names(p2)) p2[[b]] else 0
            E[a, b] <- E[a, b] + pi_k[k] * pa * pb
          }
        }
      }
    }
    E <- E * n_ok

    E_safe <- pmax(E, 1e-10)
    bvr    <- sum((as.vector(O) - as.vector(E_safe))^2 / as.vector(E_safe))
    deg_f  <- (length(cats1) - 1L) * (length(cats2) - 1L)
    p_val  <- stats::pchisq(bvr, df = max(deg_f, 1L), lower.tail = FALSE)

    rows[[idx]] <- data.frame(
      var1 = v1, var2 = v2,
      bvr = bvr, df = deg_f, p_value = p_val,
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, rows)
  result[order(-result$bvr), ]
}


# ------------------------------------------------------------------
# Internal: evaluate categorical log-density with direct effects
# ------------------------------------------------------------------
#' @keywords internal
eval_categorical_density_with_deps <- function(df, probs, dep_params,
                                               dep_pairs) {
  N  <- nrow(df)
  J  <- ncol(df)
  varnames <- colnames(df)
  ld <- numeric(N)

  # Build a lookup: child_var -> list of (parent, cond_probs_matrix) entries
  # Multiple parents for the same child are accumulated, not overwritten.
  child_map <- list()
  for (d in seq_along(dep_pairs)) {
    child_var  <- dep_pairs[[d]][2]
    parent_var <- dep_pairs[[d]][1]
    entry <- list(parent = parent_var, cond_probs = dep_params[[d]])
    if (is.null(child_map[[child_var]])) {
      child_map[[child_var]] <- list(entry)
    } else {
      child_map[[child_var]] <- c(child_map[[child_var]], list(entry))
    }
  }

  for (j in seq_len(J)) {
    vname  <- varnames[j]
    values <- as.character(df[[j]])

    if (vname %in% names(child_map)) {
      # This variable is a child; iterate over all parent effects
      dep_entries <- child_map[[vname]]
      for (entry in dep_entries) {
        parent_values <- as.character(df[[entry$parent]])
        cond_mat      <- entry$cond_probs

        for (i in seq_len(N)) {
          if (is.na(values[i])) next

          if (!is.na(parent_values[i]) &&
              parent_values[i] %in% rownames(cond_mat) &&
              values[i] %in% colnames(cond_mat)) {
            p <- cond_mat[parent_values[i], values[i]]
            ld[i] <- ld[i] + log(max(p, 1e-15))
          } else {
            # Parent missing -> fall back to marginal
            p_vec <- probs[[vname]]
            p <- if (values[i] %in% names(p_vec)) p_vec[[values[i]]] else 1e-15
            ld[i] <- ld[i] + log(max(p, 1e-15))
          }
        }
      }
    } else {
      # Standard marginal contribution
      p_vec    <- probs[[vname]]
      assigned <- p_vec[values]
      assigned[is.na(assigned)] <- NA_real_
      assigned[!is.na(assigned) & assigned < 1e-15] <- 1e-15
      contrib <- log(assigned)
      contrib[is.na(values) | is.na(contrib)] <- 0
      ld <- ld + contrib
    }
  }

  ld
}


# ------------------------------------------------------------------
# Internal: M-step for conditional probability tables
# ------------------------------------------------------------------
#' @keywords internal
update_cat_dep <- function(df, weights, dep_pairs) {
  dep_params <- vector("list", length(dep_pairs))

  for (d in seq_along(dep_pairs)) {
    parent_var <- dep_pairs[[d]][1]
    child_var  <- dep_pairs[[d]][2]

    p_vals <- as.character(df[[parent_var]])
    c_vals <- as.character(df[[child_var]])

    p_cats <- sort(unique(stats::na.omit(p_vals)))
    c_cats <- sort(unique(stats::na.omit(c_vals)))

    cond_mat <- matrix(0, length(p_cats), length(c_cats),
                       dimnames = list(p_cats, c_cats))

    for (a in p_cats) {
      p_hit  <- !is.na(p_vals) & (p_vals == a) & !is.na(c_vals)
      W_a    <- sum(weights[p_hit])

      if (W_a < 1e-15) {
        cond_mat[a, ] <- 1.0 / length(c_cats)
      } else {
        for (b in c_cats) {
          both <- p_hit & (c_vals == b)
          cond_mat[a, b] <- sum(weights[both]) / W_a
        }
      }

      # Floor and renormalise
      cond_mat[a, ] <- pmax(cond_mat[a, ], 1e-10)
      cond_mat[a, ] <- cond_mat[a, ] / sum(cond_mat[a, ])
    }

    dep_params[[d]] <- cond_mat
  }

  dep_params
}
