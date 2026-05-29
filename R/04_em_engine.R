# ==============================================================================
# File: R/04_em_engine.R
# mixLCA - core Expectation-Maximization engine.
# ==============================================================================

#' EM Engine for mixLCA
#'
#' Runs the full EM algorithm: initialises parameters via k-means,
#' iterates E/M steps until convergence, and returns a fitted
#' \code{mixLCA} object.
#'
#' @param data Data frame.
#' @param continuous Character vector of continuous indicator names, or
#'   NULL.
#' @param categorical Character vector of categorical indicator names, or
#'   NULL.
#' @param concomitant Character vector of concomitant predictor names, or
#'   NULL.
#' @param n_classes Integer >= 2.
#' @param dependence One of \code{"none"}, \code{"full"},
#'   \code{"penalized"}.
#' @param penalty Numeric L1 penalty for penalized dependence.
#' @param max_iter Maximum EM iterations.
#' @param tol Convergence tolerance on absolute log-likelihood change.
#' @param cat_direct_effects List of direct effect pairs, or NULL.
#' @param spectral_rank Integer: target rank for SLD.
#' @param spectral_pool Logical: pool Burt matrices across classes.
#' @param init_model Optional prior \code{mixLCA} object for warm-start.
#' @return An object of class \code{mixLCA}.
#' @keywords internal
run_em <- function(data, continuous = NULL, categorical = NULL,
                   concomitant = NULL, n_classes = 2L,
                   dependence = "full", penalty = 0,
                   max_iter = 500L, tol = 1e-6,
                   cat_direct_effects = NULL,
                   spectral_rank = 0L, spectral_pool = FALSE,
                   init_model = NULL, kmeans_nstart = 1L) {

  K <- as.integer(n_classes)

  has_concom <- !is.null(concomitant)
  N <- nrow(data)

  # ---- Prepare data matrices ----
  Y <- if (!is.null(continuous))
    as.matrix(data[, continuous, drop = FALSE]) else NULL
  D <- if (!is.null(categorical))
    data[, categorical, drop = FALSE] else NULL

  if (has_concom) {
    # Build the model frame once and capture its terms + factor levels.
    # Saving these on the fitted object lets predict() rebuild the
    # design matrix on newdata even when a factor level is missing in
    # the subset (which would otherwise drop a column and crash the
    # X %*% coefs multiplication).
    f_form <- if (inherits(concomitant, "formula"))
      concomitant
    else
      stats::as.formula(paste("~", paste(concomitant, collapse = " + ")))
    mf             <- stats::model.frame(f_form, data = data,
                                         na.action = stats::na.pass)
    concom_terms   <- stats::terms(mf)
    concom_xlevels <- stats::.getXlevels(concom_terms, mf)
    X              <- stats::model.matrix(concom_terms, mf)
    # Strip the .Environment attr from the saved terms: by default it
    # captures the caller's environment, which can bloat saved fits with
    # the entire calling namespace. baseenv() is safe to look up in.
    attr(concom_terms, ".Environment") <- baseenv()
  } else {
    X <- matrix(1, nrow = N, ncol = 1L)
    colnames(X) <- "(Intercept)"
    concom_terms   <- NULL
    concom_xlevels <- NULL
  }
  P <- ncol(X)

  # ---- Warm-start flag ----
  warm <- !is.null(init_model) && inherits(init_model, "mixLCA") &&
          init_model$n_classes == K

  # ---- Spectral Local Dependence Setup ----
  ranks <- as.integer(spectral_rank)
  if (length(ranks) == 1L) ranks <- rep(ranks, K)
  if (length(ranks) != K) stop("`spectral_rank` must be length 1 or K.")
  has_spectral <- !is.null(D) && any(ranks > 0L)
  if (has_spectral) {
    spec_data <- prep_spectral_data(D)
    spec_A <- list()
    spec_V <- list()
    spec_lambda <- list()
    max_d <- spec_data$C - spec_data$J
    ranks <- pmin(ranks, max_d)
    spectral_rank <- ranks  # canonical vector form

    for (k in seq_len(K)) {
      if (warm && !is.null(init_model$cat_spectral_params) &&
          length(init_model$cat_spectral_params$A_star) >= k) {
        spec_A[[k]] <- init_model$cat_spectral_params$A_star[[k]]
        # Safely pad/trim V if rank changed
        old_V <- init_model$cat_spectral_params$V_d[[k]]
        new_d <- max(ranks[k], 1L)
        new_V <- matrix(0, nrow = spec_data$C, ncol = new_d)
        if (!is.null(old_V)) {
          copy_d <- min(ncol(old_V), new_d)
          if (copy_d > 0L) new_V[, seq_len(copy_d)] <- old_V[, seq_len(copy_d)]
        }
        spec_V[[k]] <- new_V
        spec_lambda[[k]] <- init_model$cat_spectral_params$spectra[[k]]
      } else {
        spec_A[[k]] <- matrix(0, nrow = spec_data$C, ncol = spec_data$C)
        spec_V[[k]] <- matrix(0, nrow = spec_data$C, ncol = max(ranks[k], 1L))
        spec_lambda[[k]] <- rep(0, spec_data$C)
      }
    }
  }

  # ---- Initialize parameters ----
  # Random init reads from the global RNG, exactly like stats::kmeans()
  # and uwot::umap(). Users get a reproducible fit by calling
  # set.seed() before fit_lca().

  if (warm) {
    posteriors <- init_model$posteriors
    if (nrow(posteriors) != N) {
      raw        <- matrix(stats::runif(N * K), nrow = N, ncol = K)
      posteriors <- sweep(raw, 1, rowSums(raw), "/")
      warm <- FALSE
    }
  } else {
    raw        <- matrix(stats::runif(N * K), nrow = N, ncol = K)
    posteriors <- sweep(raw, 1, rowSums(raw), "/")
  }

  # Concomitant coefficients
  if (warm && !is.null(init_model$concomitant_coefs) &&
      nrow(init_model$concomitant_coefs) == P) {
    coefs <- init_model$concomitant_coefs
  } else {
    coefs <- matrix(0, nrow = P, ncol = K - 1L)
    rownames(coefs) <- colnames(X)
  }

  # Continuous parameters
  means       <- list()
  covariances <- list()
  if (!is.null(Y)) {
    if (warm && !is.null(init_model$continuous_params)) {
      means       <- init_model$continuous_params$means
      covariances <- init_model$continuous_params$covariances
    } else {
      d     <- ncol(Y)
      Y_imp <- Y
      for (col in seq_len(d)) {
        na_idx <- is.na(Y_imp[, col])
        if (any(na_idx))
          Y_imp[na_idx, col] <- mean(Y_imp[, col], na.rm = TRUE)
      }
      # kmeans can fail with empty-cluster errors on collinear or
      # heavily skewed data; fall back to a random row pick when it does.
      km <- tryCatch(
        stats::kmeans(Y_imp, centers = K, nstart = kmeans_nstart),
        error = function(e) {
          idx <- sample.int(nrow(Y_imp), K, replace = FALSE)
          list(centers = Y_imp[idx, , drop = FALSE])
        }
      )
      for (k in seq_len(K)) {
        means[[k]]       <- km$centers[k, ]
        v_vec            <- apply(Y_imp, 2, stats::var, na.rm = TRUE)
        # diag(v) with length(v)==1 allocates a v x v identity, hence
        # explicit nrow.
        covariances[[k]] <- diag(v_vec, nrow = length(v_vec))
      }
    }
  }

  # Categorical parameters
  cat_params <- list()
  if (!is.null(D)) {
    if (warm && !is.null(init_model$categorical_params)) {
      cat_params <- init_model$categorical_params
    } else {
      for (k in seq_len(K))
        cat_params[[k]] <- update_categorical(D, posteriors[, k])
    }
  }

  # Direct effect conditional probability tables
  has_deps <- !is.null(cat_direct_effects) && length(cat_direct_effects) > 0
  cat_dep_params <- list()
  if (has_deps && !is.null(D)) {
    for (k in seq_len(K))
      cat_dep_params[[k]] <- update_cat_dep(D, posteriors[, k],
                                            cat_direct_effects)
  }

  # ---- EM iterations ----
  ll_history <- numeric(max_iter)
  converged  <- FALSE
  ll         <- -Inf
  # Declare in this lexical scope so the EM loop never has to rely on
  # exists(), which would search the entire call stack.
  active_offdiag_list <- NULL

  for (iter in seq_len(max_iter)) {

    # === E-step ===

    priors    <- compute_priors(X, coefs)
    log_dens  <- matrix(0, nrow = N, ncol = K)

    for (k in seq_len(K)) {
      if (!is.null(Y))
        log_dens[, k] <- log_dens[, k] +
          eval_continuous_density(Y, means[[k]], covariances[[k]])
      if (!is.null(D)) {
        if (has_spectral && ranks[k] > 0L) {
          pi_c <- flatten_cat_params(cat_params[[k]], spec_data$item_indices, spec_data$C)
          log_dens[, k] <- log_dens[, k] +
            eval_spectral_density(spec_data, pi_c, spec_A[[k]])
        } else if (has_spectral && ranks[k] == 0L) {
          # Class k has no SLD: use plain marginal density
          log_dens[, k] <- log_dens[, k] +
            eval_categorical_density(D, cat_params[[k]])
        } else if (has_deps) {
          log_dens[, k] <- log_dens[, k] +
            eval_categorical_density_with_deps(
              D, cat_params[[k]], cat_dep_params[[k]],
              cat_direct_effects)
        } else {
          log_dens[, k] <- log_dens[, k] +
            eval_categorical_density(D, cat_params[[k]])
        }
      }
    }

    log_joint <- log(priors + 1e-300) + log_dens
    log_marg  <- apply(log_joint, 1, log_sum_exp)

    # Guard against NaN/Inf from degenerate classes
    bad <- !is.finite(log_marg)
    if (any(bad)) log_marg[bad] <- -1e10

    posteriors <- exp(log_joint - log_marg)
    ll         <- sum(log_marg)
    ll_history[iter] <- ll

    if (K > 1L) {
      alpha_prior <- 1.01
      posteriors  <- (posteriors * N + (alpha_prior - 1)) /
                       (N + K * (alpha_prior - 1))
      posteriors  <- sweep(posteriors, 1, rowSums(posteriors), "/")
    }

    # Convergence check
    if (iter > 1L &&
        abs(ll - ll_history[iter - 1L]) < tol) {
      converged  <- TRUE
      ll_history <- ll_history[seq_len(iter)]
      break
    }

    # === M-step ===

    if (has_spectral) {
      spec_A_old <- spec_A
      Sigma_list <- list()
    }

    for (k in seq_len(K)) {
      wk <- posteriors[, k]

      if (!is.null(Y)) {
        upd             <- update_continuous(Y, wk, means[[k]],
                                            covariances[[k]],
                                            dependence = dependence,
                                            penalty    = penalty)
        means[[k]]       <- upd$mean
        covariances[[k]] <- upd$covariance
        if (!is.null(upd$active_offdiag)) {
          if (is.null(active_offdiag_list)) active_offdiag_list <- vector("list", K)
          active_offdiag_list[[k]] <- upd$active_offdiag
        }
      }

      if (!is.null(D)) {
        cat_params[[k]] <- update_categorical(D, wk)

        if (has_spectral) {
          pi_c <- flatten_cat_params(cat_params[[k]], spec_data$item_indices, spec_data$C)
          E_c <- matrix(pi_c, nrow = N, ncol = spec_data$C, byrow = TRUE)
          Z_imp <- spec_data$Z
          Z_imp[spec_data$Z_mis] <- E_c[spec_data$Z_mis]
          R_c <- Z_imp - E_c
          W_sum <- sum(wk)
          if (W_sum > 1e-15) {
            R_w <- sweep(R_c, 1, sqrt(wk), "*")
            Sigma_list[[k]] <- crossprod(R_w) / W_sum
          } else {
            Sigma_list[[k]] <- matrix(0, spec_data$C, spec_data$C)
          }
        }
      }

      if (has_deps && !is.null(D) && !has_spectral)
        cat_dep_params[[k]] <- update_cat_dep(D, wk, cat_direct_effects)
    }

    if (K > 1L) {
      coefs <- update_concomitant(X, posteriors, coefs)
    }

    if (has_spectral) {
      if (isTRUE(spectral_pool)) {
        Sigma_pool <- matrix(0, nrow = spec_data$C, ncol = spec_data$C)
        for (k in seq_len(K)) Sigma_pool <- Sigma_pool + (sum(posteriors[, k]) / N) * Sigma_list[[k]]
        res <- spectral_projector(Sigma_pool, max(spectral_rank), spec_data$M)
        for (k in seq_len(K)) { spec_A[[k]] <- res$A_star; spec_V[[k]] <- res$V_d; spec_lambda[[k]] <- res$lambda }
      } else {
        for (k in seq_len(K)) {
          if (spectral_rank[k] > 0L) {
            res <- spectral_projector(Sigma_list[[k]], spectral_rank[k], spec_data$M)
            spec_A[[k]] <- res$A_star; spec_V[[k]] <- res$V_d; spec_lambda[[k]] <- res$lambda
          } else {
            # d=0: no spectral shift for this class
            spec_A[[k]] <- matrix(0, nrow = spec_data$C, ncol = spec_data$C)
          }
        }
      }

      # Generalized EM Step-halving safeguard
      halves <- 0L
      while (halves < 5L) {
        test_log_dens <- matrix(0, nrow = N, ncol = K)
        test_priors <- compute_priors(X, coefs)
        for (k in seq_len(K)) {
          if (!is.null(Y)) test_log_dens[, k] <- test_log_dens[, k] + eval_continuous_density(Y, means[[k]], covariances[[k]])
          if (!is.null(D)) {
            if (spectral_rank[k] > 0L) {
              pi_c <- flatten_cat_params(cat_params[[k]], spec_data$item_indices, spec_data$C)
              test_log_dens[, k] <- test_log_dens[, k] + eval_spectral_density(spec_data, pi_c, spec_A[[k]])
            } else {
              test_log_dens[, k] <- test_log_dens[, k] + eval_categorical_density(D, cat_params[[k]])
            }
          }
        }
        test_log_joint <- log(test_priors + 1e-300) + test_log_dens
        test_ll <- sum(apply(test_log_joint, 1, log_sum_exp))
        if (test_ll >= ll - 1e-6) break else {
          for (k in seq_len(K)) spec_A[[k]] <- 0.5 * (spec_A[[k]] + spec_A_old[[k]])
          halves <- halves + 1L
        }
      }

      # If step-halving occurred, re-decompose to synchronize V and lambda
      if (halves > 0L) {
        for (k in seq_len(K)) {
          d_k <- if (isTRUE(spectral_pool)) max(spectral_rank) else spectral_rank[k]
          if (d_k > 0L) {
            eig_h <- eigen(spec_A[[k]] * spec_data$M, symmetric = TRUE)
            spec_V[[k]] <- eig_h$vectors[, seq_len(d_k), drop = FALSE]
            spec_lambda[[k]] <- eig_h$values
          }
        }
      }
    }
  }

  if (!converged)
    ll_history <- ll_history[seq_len(iter)]

  # ---- Assemble mixLCA object ----
  result <- list(
    n_classes          = K,
    log_lik            = ll,
    posteriors         = posteriors,
    continuous_params  = if (!is.null(Y))
      list(means = means, covariances = covariances,
           active_offdiag = active_offdiag_list) else NULL,
    categorical_params = if (!is.null(D)) cat_params else NULL,
    cat_dep_params     = if (has_deps && !has_spectral) cat_dep_params else NULL,
    cat_spectral_params = if (has_spectral) list(A_star = spec_A, V_d = spec_V, spectra = spec_lambda, encoding = spec_data) else NULL,
    concomitant_coefs  = if (has_concom) coefs else NULL,
    n_obs              = N,
    n_obs_effective    = NA_integer_,
    convergence        = list(
      iterations = length(ll_history),
      converged  = converged,
      ll_history = ll_history
    ),
    specs = list(
      continuous         = continuous,
      categorical        = categorical,
      concomitant        = concomitant,
      dependence         = dependence,
      penalty            = penalty,
      cat_direct_effects = cat_direct_effects,
      spectral_rank      = if (has_spectral) as.integer(spectral_rank) else rep(0L, K),
      spectral_pool      = if (has_spectral) isTRUE(spectral_pool) else FALSE,
      concom_terms       = concom_terms,
      concom_xlevels     = concom_xlevels
    )
  )
  class(result) <- "mixLCA"

  # Effective N: rows contributing at least one non-missing indicator
  indicator_cols <- c(continuous, categorical)
  indicator_df <- data[, indicator_cols, drop = FALSE]
  n_eff <- sum(rowSums(!is.na(indicator_df)) > 0L)
  result$n_obs_effective <- n_eff

  result$n_params <- count_params(result)
  result$loglik   <- result$log_lik  # alias for convenience

  # Pre-compute and cache fit indices on the object
  fi <- fit_indices(result)
  result$AIC     <- fi$AIC
  result$BIC     <- fi$BIC
  result$aBIC    <- fi$aBIC
  result$ICL     <- fi$ICL
  result$entropy <- fi$entropy

  result
}
