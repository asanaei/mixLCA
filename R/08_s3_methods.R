# ==============================================================================
# File: R/08_s3_methods.R
# mixLCA - standard S3 methods: print, summary, predict, coef, logLik, nobs
# ==============================================================================

#' Extract Posterior Class Probabilities
#'
#' Accessor returning the N x K posterior probability matrix. Use this
#' rather than \code{model$posteriors} so downstream code is insulated
#' against future restructuring of internal fields.
#'
#' @param x A fitted \code{mixLCA} object.
#' @param ... Unused.
#' @return Numeric matrix.
#'
#' @examples
#' \donttest{
#' data(voter_perceptions)
#' fit <- fit_lca(voter_perceptions,
#'                categorical = names(voter_perceptions),
#'                n_classes   = 2,
#'                control     = lca_control(n_starts = 2, seed = 110),
#'                verbose     = FALSE)
#' post <- get_posteriors(fit)
#' dim(post)
#' head(round(post, 3))
#' colSums(post) / nrow(post)   # estimated class proportions
#' }
#' @export
get_posteriors <- function(x, ...) UseMethod("get_posteriors")

#' @rdname get_posteriors
#' @export
get_posteriors.mixLCA <- function(x, ...) x$posteriors

#' Extract Spectral Loadings
#'
#' Accessor returning the SLD loadings table (or NULL when the model
#' does not use Spectral Local Dependence). Use this rather than
#' \code{model$cat_spectral_params} so downstream code is insulated
#' against future restructuring of internal fields.
#'
#' @param x A fitted \code{mixLCA} object.
#' @param ... Unused.
#' @return A data frame or NULL.
#'
#' @examples
#' \donttest{
#' data(voter_perceptions)
#' fit <- fit_lca(voter_perceptions,
#'                categorical   = names(voter_perceptions),
#'                n_classes     = 2,
#'                spectral_rank = c(1L, 1L),
#'                control       = lca_control(n_starts = 2, seed = 110),
#'                verbose       = FALSE)
#' ld <- get_loadings(fit)
#' head(ld, 12)
#' }
#' @export
get_loadings <- function(x, ...) UseMethod("get_loadings")

#' @rdname get_loadings
#' @export
get_loadings.mixLCA <- function(x, ...) {
  if (!is.null(x$cat_spectral_params)) spectral_loadings(x) else NULL
}

#' Print a mixLCA Object
#'
#' @param x A \code{mixLCA} object.
#' @param ... Unused.
#'
#' @examples
#' \donttest{
#' data(voter_perceptions)
#' fit <- fit_lca(voter_perceptions,
#'                categorical = names(voter_perceptions),
#'                n_classes   = 2,
#'                control     = lca_control(n_starts = 2, seed = 110),
#'                verbose     = FALSE)
#' print(fit)
#' }
#' @export
print.mixLCA <- function(x, ...) {
  cat("\nLatent Class Model - mixLCA\n")
  cat("===========================\n")
  cat("Classes        :", x$n_classes, "\n")
  is_mixed <- !is.null(x$specs$continuous) && isTRUE(any(x$specs$spectral_rank > 0L))
  is_comp  <- is.null(x$specs$continuous) && isTRUE(any(x$specs$spectral_rank > 0L))
  ll_label <- if (is_mixed) "Partial Composite LL " else if (is_comp) "Composite Log-Lik.   " else "Log-likelihood       "
  cat(sprintf("%s:", ll_label), round(x$log_lik, 4), "\n")
  cat("Parameters     :", x$n_params, "\n")
  n_eff <- if (!is.null(x$n_obs_effective) && !is.na(x$n_obs_effective))
             x$n_obs_effective else x$n_obs
  cat("Observations   :", n_eff, "\n")
  if (n_eff != x$n_obs)
    cat("  (total rows  :", x$n_obs, ", rows with missing indicators excluded from N)\n")
  cat("Converged      :", x$convergence$converged,
      paste0("(", x$convergence$iterations, " iterations)"), "\n")

  cat("\nSpecification:\n")
  if (!is.null(x$specs$continuous))
    cat("  Continuous indicators :", length(x$specs$continuous), "\n")
  if (!is.null(x$specs$categorical))
    cat("  Categorical indicators:", length(x$specs$categorical), "\n")
  if (!is.null(x$specs$concomitant))
    cat("  Concomitant predictors:", length(x$specs$concomitant), "\n")
  cat("  Dependence structure  :", x$specs$dependence, "\n")

  if (isTRUE(any(x$specs$spectral_rank > 0L))) {
    cat("  Categorical dep.      : Spectral Local Dependence (rank d = ", paste(x$specs$spectral_rank, collapse = ", "), ")\n", sep = "")
    cat("  Spectral basis        :", if (isTRUE(x$specs$spectral_pool)) "pooled" else "class-specific", "\n")
  }

  cat("\nClass proportions:\n")
  props        <- colMeans(x$posteriors)
  names(props) <- paste("Class", seq_len(x$n_classes))
  print(round(props, 4))
  cat("\n")
  invisible(x)
}

#' Summary Method for mixLCA
#'
#' Prints fit indices, class proportions, measurement parameters, and
#' concomitant coefficients with standard errors when \code{data} is
#' supplied.
#'
#' @param object A \code{mixLCA} object.
#' @param data Optional data frame: required for concomitant standard
#'   errors.
#' @param ... Unused.
#'
#' @examples
#' \donttest{
#' data(health_screening)
#' fit <- fit_lca(health_screening,
#'                continuous  = c("marker_1","marker_2","marker_3","marker_4"),
#'                concomitant = ~ age,
#'                n_classes   = 2,
#'                control     = lca_control(n_starts = 3, seed = 110),
#'                verbose     = FALSE)
#' summary(fit, data = health_screening)
#' }
#' @export
summary.mixLCA <- function(object, data = NULL, ...) {
  fi <- fit_indices(object)
  K  <- object$n_classes

  cat("\nSummary - mixLCA\n")
  cat("================\n\n")

  cat("Fit Indices:\n")
  is_mixed <- !is.null(object$specs$continuous) && isTRUE(any(object$specs$spectral_rank > 0L))
  is_comp  <- is.null(object$specs$continuous) && isTRUE(any(object$specs$spectral_rank > 0L))
  ll_label <- if (is_mixed) "Partial Composite LL" else if (is_comp) "Composite Log-Lik. " else "Log-likelihood     "
  cat(sprintf("  %s:", ll_label), round(fi$log_lik, 4), "\n")
  bic_label <- if (isTRUE(fi$is_composite)) "  BIC (naive)       :" else "  BIC               :"
  cat("  AIC               :", round(fi$AIC, 4), "\n")
  cat(bic_label, round(fi$BIC, 4), "\n")
  cat("  Sample-adj. BIC   :", round(fi$aBIC, 4), "\n")
  cat("  Entropy           :", round(fi$entropy, 4), "\n")
  cat("  ICL               :", round(fi$ICL, 4), "\n")
  if (isTRUE(fi$is_composite))
    cat("  Note: BIC is naive under composite likelihood; prefer V-fold CV for rank selection.\n")
  cat("\n")

  cat("Class Proportions:\n")
  props        <- colMeans(object$posteriors)
  names(props) <- paste("Class", seq_len(K))
  print(round(props, 4))
  cat("\n")

  # Continuous measurement parameters
  if (!is.null(object$continuous_params)) {
    cont <- object$specs$continuous
    cat("Continuous Indicator Means:\n")
    mu_mat           <- do.call(rbind, object$continuous_params$means)
    rownames(mu_mat) <- paste("Class", seq_len(K))
    colnames(mu_mat) <- cont
    print(round(mu_mat, 4))
    cat("\n")

    if (object$specs$dependence != "none") {
      for (k in seq_len(K)) {
        cat("Covariance (Class ", k, "):\n", sep = "")
        Sk             <- object$continuous_params$covariances[[k]]
        rownames(Sk)   <- cont
        colnames(Sk)   <- cont
        print(round(Sk, 4))
        cat("\n")
      }
    }
  }

  # Categorical measurement parameters
  if (!is.null(object$categorical_params)) {
    cat("Categorical Response Probabilities:\n")
    for (k in seq_len(K)) {
      cat("  Class ", k, ":\n", sep = "")
      for (vname in names(object$categorical_params[[k]])) {
        probs <- object$categorical_params[[k]][[vname]]
        cat("    ", vname, ": ",
            paste(names(probs), "=", round(probs, 3), collapse = "  "),
            "\n")
      }
    }
    cat("\n")
  }

  # Direct effects (categorical local dependence)
  if (!is.null(object$cat_dep_params) &&
      !is.null(object$specs$cat_direct_effects)) {
    deps <- object$specs$cat_direct_effects
    cat("Categorical Direct Effects (local dependence):\n")
    for (d in seq_along(deps)) {
      parent <- deps[[d]][1]
      child  <- deps[[d]][2]
      cat("  ", parent, " -> ", child, ":\n", sep = "")
      for (k in seq_len(K)) {
        cond_mat <- object$cat_dep_params[[k]][[d]]
        cat("    Class ", k, ":\n", sep = "")
        cat("    P(", child, " | ", parent, "):\n", sep = "")
        print_mat <- round(cond_mat, 3)
        for (rn in rownames(print_mat)) {
          cat("      ", parent, "=", rn, ": ",
              paste(colnames(print_mat), "=", print_mat[rn, ],
                    collapse = "  "), "\n")
        }
      }
    }
    cat("\n")
  }

  if (isTRUE(any(object$specs$spectral_rank > 0L)) && !is.null(object$cat_spectral_params)) {
    sp <- object$cat_spectral_params
    cat("Spectral Basis:\n")
    d_vec <- object$specs$spectral_rank
    if (length(d_vec) == 1L) d_vec <- rep(d_vec, K)

    if (isTRUE(object$specs$spectral_pool)) {
      eigs <- sp$spectra[[1L]]; n_show <- min(max(d_vec) + 2L, sum(eigs > 1e-10))
      cat("  Pooled eigenvalues (first ", n_show, "): ", paste(round(utils::head(eigs, n_show), 4), collapse = ", "), "\n\n", sep = "")
    } else {
      for (k in seq_len(K)) {
        if (d_vec[k] == 0L) {
          cat("  Class ", k, ": d = 0 (local independence)\n", sep = "")
        } else {
          eigs <- sp$spectra[[k]]; n_show <- min(d_vec[k] + 2L, sum(eigs > 1e-10))
          cat("  Class ", k, " eigenvalues (first ", n_show, "): ", paste(round(utils::head(eigs, n_show), 4), collapse = ", "), "\n", sep = "")
        }
      }
      cat("\n")
    }
  }

  # Concomitant coefficients
  if (!is.null(object$concomitant_coefs)) {
    cat("Concomitant Coefficients (reference = Class 1):\n")
    pred_names <- rownames(object$concomitant_coefs)

    se_mat <- NULL
    if (!is.null(data))
      se_mat <- concomitant_se(object, data)

    for (k in seq_len(K - 1L)) {
      cat("\n  -> Class ", k + 1L, ":\n", sep = "")
      est <- object$concomitant_coefs[, k]

      if (!is.null(se_mat)) {
        se  <- se_mat[, k]
        z   <- est / se
        pv  <- 2 * stats::pnorm(-abs(z))
        tab <- cbind(Estimate = est, SE = se, z = z, p = pv)
        rownames(tab) <- pred_names
        stats::printCoefmat(tab, digits = 4, signif.stars = TRUE,
                            P.values = TRUE, has.Pvalue = TRUE)
      } else {
        names(est) <- pred_names
        print(round(est, 4))
      }
    }
    cat("\n")
    if (is.null(data) && !is.null(object$specs$concomitant))
      cat("  (Pass `data` to summary() for standard errors.)\n\n")
  }

  invisible(object)
}

#' Predict Method for mixLCA
#'
#' Returns posterior class probabilities, modal class assignments, or
#' the full diagnostic dataframe for each observation. When
#' \code{newdata} is supplied, posteriors are computed from scratch
#' using the estimated model parameters, including Spectral Local
#' Dependence shifts if applicable.
#'
#' Rows of \code{newdata} with missing concomitant values are not
#' scored and are returned as \code{NA} so that the output length
#' always equals \code{nrow(newdata)} (the standard \code{na.exclude}
#' contract). This makes \code{cbind(newdata, predict(model, newdata))}
#' safe.
#'
#' @param object A \code{mixLCA} object.
#' @param newdata Optional data frame for out-of-sample prediction.
#' @param type One of \code{"prob"} (default; N x K matrix of posterior
#'   probabilities), \code{"class"} (integer vector of modal class
#'   assignments), or \code{"all"} (legacy data frame containing
#'   \code{P_class_*}, \code{modal_class}, \code{max_posterior},
#'   and -- for out-of-sample -- \code{log_lik}).
#' @param ... Unused.
#' @return A matrix, integer vector, or data frame, per \code{type}.
#'
#' @examples
#' \donttest{
#' data(voter_perceptions)
#' fit <- fit_lca(voter_perceptions,
#'                categorical = names(voter_perceptions),
#'                n_classes   = 3,
#'                control     = lca_control(n_starts = 2, seed = 110),
#'                verbose     = FALSE)
#'
#' # Default: N x K posterior probability matrix
#' head(predict(fit))
#'
#' # Modal class assignments as an integer vector
#' head(predict(fit, type = "class"))
#'
#' # Full diagnostic data frame
#' head(predict(fit, type = "all"))
#'
#' # Out-of-sample prediction; rows with NA in concomitants pad with NA
#' new_rows <- voter_perceptions[1:6, ]
#' predict(fit, newdata = new_rows, type = "class")
#' }
#' @export
predict.mixLCA <- function(object, newdata = NULL,
                           type = c("prob", "class", "all"), ...) {

  type <- match.arg(type)
  K    <- object$n_classes

  if (is.null(newdata)) {
    post <- object$posteriors
    if (is.null(colnames(post)))
      colnames(post) <- paste0("P_class_", seq_len(K))
    if (type == "prob")  return(post)
    if (type == "class") return(as.integer(apply(post, 1, which.max)))
    out              <- as.data.frame(post)
    colnames(out)    <- paste0("P_class_", seq_len(K))
    out$modal_class  <- apply(post, 1, which.max)
    out$max_posterior <- apply(post, 1, max)
    return(out)
  }

  N_in   <- nrow(newdata)
  cont   <- object$specs$continuous
  cat_   <- object$specs$categorical
  concom <- object$specs$concomitant

  # Determine which rows are scorable (concomitant predictors complete);
  # rows with NA in concomitant inputs are returned as NA.
  if (!is.null(concom)) {
    concom_vars <- if (inherits(concom, "formula")) all.vars(concom) else concom
    scorable    <- stats::complete.cases(newdata[, concom_vars, drop = FALSE])
  } else {
    scorable <- rep(TRUE, N_in)
  }

  na_rows <- which(!scorable)
  in_data <- newdata[scorable, , drop = FALSE]
  N       <- nrow(in_data)

  if (N == 0L) {
    # All rows missing concomitant; return all-NA output of correct shape.
    post_full <- matrix(NA_real_, nrow = N_in, ncol = K)
    colnames(post_full) <- paste0("P_class_", seq_len(K))
    if (type == "prob")  return(post_full)
    if (type == "class") return(rep(NA_integer_, N_in))
    out_full <- as.data.frame(post_full)
    out_full$modal_class   <- NA_integer_
    out_full$max_posterior <- NA_real_
    out_full$log_lik       <- NA_real_
    return(out_full)
  }

  Y <- if (!is.null(cont)) as.matrix(in_data[, cont, drop = FALSE]) else NULL
  D <- if (!is.null(cat_)) in_data[, cat_, drop = FALSE] else NULL

  if (!is.null(concom) && !is.null(object$concomitant_coefs)) {
    if (inherits(concom, "formula")) {
      X <- stats::model.matrix(concom, data = in_data)
    } else {
      f_str <- paste("~", paste(concom, collapse = " + "))
      X <- stats::model.matrix(stats::as.formula(f_str), data = in_data)
    }
    priors <- compute_priors(X, object$concomitant_coefs)
  } else {
    coefs_dummy <- matrix(0, nrow = 1L, ncol = K - 1L)
    rownames(coefs_dummy) <- "(Intercept)"
    X <- matrix(1, nrow = N, ncol = 1L)
    priors <- compute_priors(X, coefs_dummy)
  }

  has_spectral <- !is.null(D) && isTRUE(any(object$specs$spectral_rank > 0L)) &&
                  !is.null(object$cat_spectral_params)
  if (has_spectral) {
    spec_enc      <- object$cat_spectral_params$encoding
    spec_data_new <- encode_newdata_spectral(D, spec_enc)
  }

  has_deps <- !is.null(object$cat_dep_params) &&
              !is.null(object$specs$cat_direct_effects)

  log_dens <- matrix(0, nrow = N, ncol = K)
  for (k in seq_len(K)) {
    if (!is.null(Y))
      log_dens[, k] <- log_dens[, k] +
        eval_continuous_density(Y, object$continuous_params$means[[k]],
                                object$continuous_params$covariances[[k]])
    if (!is.null(D)) {
      d_vec <- object$specs$spectral_rank
      if (length(d_vec) == 1L) d_vec <- rep(d_vec, K)
      if (has_spectral && d_vec[k] > 0L) {
        pi_c <- flatten_cat_params(object$categorical_params[[k]],
                                   spec_enc$item_indices, spec_enc$C)
        log_dens[, k] <- log_dens[, k] +
          eval_spectral_density_oos(spec_data_new, spec_enc, pi_c,
                                    object$cat_spectral_params$A_star[[k]])
      } else if (has_deps) {
        log_dens[, k] <- log_dens[, k] +
          eval_categorical_density_with_deps(
            D, object$categorical_params[[k]],
            object$cat_dep_params[[k]],
            object$specs$cat_direct_effects)
      } else {
        log_dens[, k] <- log_dens[, k] +
          eval_categorical_density(D, object$categorical_params[[k]])
      }
    }
  }

  log_joint <- log(priors + 1e-300) + log_dens
  log_marg  <- apply(log_joint, 1, log_sum_exp)
  post      <- exp(log_joint - log_marg)

  # Pad back to nrow(newdata)
  post_full <- matrix(NA_real_, nrow = N_in, ncol = K)
  post_full[scorable, ] <- post
  colnames(post_full)   <- paste0("P_class_", seq_len(K))

  if (type == "prob") return(post_full)

  modal_full <- rep(NA_integer_, N_in)
  modal_full[scorable] <- as.integer(apply(post, 1, which.max))

  if (type == "class") return(modal_full)

  max_full <- rep(NA_real_, N_in)
  max_full[scorable] <- apply(post, 1, max)
  ll_full  <- rep(NA_real_, N_in)
  ll_full[scorable]  <- log_marg

  out_full <- as.data.frame(post_full)
  out_full$modal_class   <- modal_full
  out_full$max_posterior <- max_full
  out_full$log_lik       <- ll_full
  out_full
}

#' Coef Method for mixLCA
#'
#' Returns a list of all estimated parameters.
#'
#' @param object A \code{mixLCA} object.
#' @param ... Unused.
#' @return Named list.
#'
#' @examples
#' \donttest{
#' data(health_screening)
#' fit <- fit_lca(health_screening,
#'                continuous  = c("marker_1","marker_2","marker_3","marker_4"),
#'                concomitant = ~ age,
#'                n_classes   = 2,
#'                control     = lca_control(n_starts = 2, seed = 110),
#'                verbose     = FALSE)
#' cc <- coef(fit)
#' names(cc)
#' cc$concomitant
#' cc$means
#' }
#' @export
coef.mixLCA <- function(object, ...) {
  parts <- list()

  if (!is.null(object$concomitant_coefs)) {
    cc             <- object$concomitant_coefs
    colnames(cc)   <- paste0("Class_", 2:object$n_classes)
    parts$concomitant <- cc
  }

  if (!is.null(object$continuous_params)) {
    parts$means       <- object$continuous_params$means
    parts$covariances <- object$continuous_params$covariances
  }

  if (!is.null(object$categorical_params))
    parts$categorical <- object$categorical_params

  parts
}

#' logLik Method for mixLCA
#'
#' @param object A \code{mixLCA} object.
#' @param ... Unused.
#' @return Object of class \code{logLik}.
#'
#' @examples
#' \donttest{
#' data(voter_perceptions)
#' fit <- fit_lca(voter_perceptions,
#'                categorical = names(voter_perceptions),
#'                n_classes   = 2,
#'                control     = lca_control(n_starts = 2, seed = 110),
#'                verbose     = FALSE)
#' logLik(fit)
#' }
#' @export
logLik.mixLCA <- function(object, ...) {
  ll            <- object$log_lik
  attr(ll, "df")   <- object$n_params
  n_eff <- if (!is.null(object$n_obs_effective) && !is.na(object$n_obs_effective))
             object$n_obs_effective else object$n_obs
  attr(ll, "nobs") <- n_eff
  class(ll) <- "logLik"
  ll
}

#' nobs Method for mixLCA
#'
#' @param object A \code{mixLCA} object.
#' @param ... Unused.
#' @return Integer.
#'
#' @examples
#' \donttest{
#' data(voter_perceptions)
#' fit <- fit_lca(voter_perceptions,
#'                categorical = names(voter_perceptions),
#'                n_classes   = 2,
#'                control     = lca_control(n_starts = 2, seed = 110),
#'                verbose     = FALSE)
#' nobs(fit)
#' }
#' @export
nobs.mixLCA <- function(object, ...) {
  if (!is.null(object$n_obs_effective) && !is.na(object$n_obs_effective))
    object$n_obs_effective else object$n_obs
}

#' AIC Method for mixLCA
#'
#' @param object A \code{mixLCA} object.
#' @param ... Unused.
#' @param k Numeric penalty per parameter (default 2).
#' @return Numeric AIC value.
#'
#' @examples
#' \donttest{
#' data(voter_perceptions)
#' fit <- fit_lca(voter_perceptions,
#'                categorical = names(voter_perceptions),
#'                n_classes   = 2,
#'                control     = lca_control(n_starts = 2, seed = 110),
#'                verbose     = FALSE)
#' AIC(fit)
#' # BIC-like penalty (k = log(N))
#' AIC(fit, k = log(nobs(fit)))
#' }
#' @export
AIC.mixLCA <- function(object, ..., k = 2) {
  if (!is.null(object$AIC) && k == 2) return(object$AIC)
  -2 * object$log_lik + k * object$n_params
}

#' BIC Method for mixLCA
#'
#' @param object A \code{mixLCA} object.
#' @param ... Unused.
#' @return Numeric BIC value.
#'
#' @examples
#' \donttest{
#' data(voter_perceptions)
#' fit <- fit_lca(voter_perceptions,
#'                categorical = names(voter_perceptions),
#'                n_classes   = 2,
#'                control     = lca_control(n_starts = 2, seed = 110),
#'                verbose     = FALSE)
#' BIC(fit)
#' }
#' @export
BIC.mixLCA <- function(object, ...) {
  if (!is.null(object$BIC)) return(object$BIC)
  N <- if (!is.null(object$n_obs_effective) && !is.na(object$n_obs_effective))
         object$n_obs_effective else object$n_obs
  -2 * object$log_lik + object$n_params * log(N)
}
