# ==============================================================================
# File: R/14_spectral_dependence.R
# mixLCA - Spectral Local Dependence for Categorical Indicators
# ==============================================================================

#' Prepare Categorical Data for Spectral Local Dependence
#' @keywords internal
prep_spectral_data <- function(df) {
  J <- ncol(df)
  N <- nrow(df)
  varnames <- colnames(df)
  
  cat_levels <- list()
  C <- 0L
  item_indices <- list()
  
  for (v in varnames) {
    cats <- sort(unique(stats::na.omit(as.character(df[[v]]))))
    if (length(cats) < 2L) stop(sprintf("Categorical indicator '%s' has fewer than 2 levels.", v))
    cat_levels[[v]] <- cats
    K_j <- length(cats)
    item_indices[[v]] <- seq(C + 1L, C + K_j)
    C <- C + K_j
  }
  
  Z <- matrix(0, nrow = N, ncol = C)
  Z_mis <- matrix(FALSE, nrow = N, ncol = C)
  col_names <- character(C)
  item_of <- integer(C)
  
  for (j in seq_len(J)) {
    v <- varnames[j]
    vals <- as.character(df[[v]])
    cats <- cat_levels[[v]]
    idx <- item_indices[[v]]
    
    mis <- is.na(vals)
    for (k in seq_along(cats)) {
      Z[!mis, idx[k]] <- as.numeric(vals[!mis] == cats[k])
      Z_mis[mis, idx[k]] <- TRUE
      col_names[idx[k]] <- paste(v, cats[k], sep = "::")
      item_of[idx[k]] <- j
    }
  }
  colnames(Z) <- col_names
  
  M <- matrix(1, nrow = C, ncol = C)
  for (idx in item_indices) M[idx, idx] <- 0
  
  list(Z = Z, Z_mis = Z_mis, C = C, J = J, 
       item_indices = item_indices, cat_levels = cat_levels, M = M,
       item_names = varnames, item_of = item_of)
}

#' Flatten Categorical Parameters List
#' @keywords internal
flatten_cat_params <- function(cat_params_k, item_indices, C) {
  pi_c <- numeric(C)
  for (vname in names(item_indices)) {
    pi_c[item_indices[[vname]]] <- cat_params_k[[vname]]
  }
  pi_c
}

#' Evaluate Spectral Composite Log-Density
#' @keywords internal
eval_spectral_density <- function(spec_data, pi_c, A_star) {
  # Derive item layout vectors for the C++ engine (0-based)
  items <- names(spec_data$item_indices)
  item_starts <- vapply(items, function(v) spec_data$item_indices[[v]][1L] - 1L,
                        integer(1))
  item_sizes  <- vapply(items, function(v) length(spec_data$item_indices[[v]]),
                        integer(1))

  as.vector(eval_spectral_density_cpp(
    spec_data$Z,
    spec_data$Z_mis * 1L,     # logical -> integer for Rcpp umat
    as.numeric(pi_c),
    A_star,
    as.integer(item_starts),
    as.integer(item_sizes)
  ))
}

#' Spectral Projector Calculation
#' @keywords internal
spectral_projector <- function(Sigma_c, d, M) {
  C <- ncol(Sigma_c)
  V_d <- matrix(0, C, max(d, 1L))
  
  if (d < 1) return(list(A_star = matrix(0, C, C), V_d = V_d, lambda = rep(0, C)))
  
  Sigma_c <- (Sigma_c + t(Sigma_c)) / 2
  eig <- eigen(Sigma_c, symmetric = TRUE)
  
  d_eff <- min(d, sum(eig$values > 1e-10))
  if (d_eff > 0) {
    V_d[, seq_len(d_eff)] <- eig$vectors[, seq_len(d_eff)]
    A_c <- tcrossprod(V_d[, seq_len(d_eff), drop = FALSE])
  } else {
    A_c <- matrix(0, C, C)
  }
  
  list(A_star = A_c * M, V_d = V_d, lambda = eig$values)
}

#' Re-Encode New Data Using Training Spectral Encoding
#'
#' Builds the indicator matrix Z and missingness mask for out-of-sample
#' data, using the category levels fixed at training time. Unseen
#' categories are treated as missing.
#'
#' @param df_new Data frame of categorical indicators.
#' @param encoding Training-time encoding list from
#'   \code{model$cat_spectral_params$encoding}.
#' @return List with elements \code{Z} and \code{Z_mis}.
#' @keywords internal
encode_newdata_spectral <- function(df_new, encoding) {
  N <- nrow(df_new)
  C <- encoding$C
  Z <- matrix(0, nrow = N, ncol = C)
  Z_mis <- matrix(FALSE, nrow = N, ncol = C)

  for (v in encoding$item_names) {
    vals <- as.character(df_new[[v]])
    cats <- encoding$cat_levels[[v]]
    idx  <- encoding$item_indices[[v]]

    for (i in seq_len(N)) {
      if (is.na(vals[i]) || !(vals[i] %in% cats)) {
        Z_mis[i, idx] <- TRUE
      } else {
        hit <- which(cats == vals[i])
        Z[i, idx[hit]] <- 1
      }
    }
  }

  list(Z = Z, Z_mis = Z_mis)
}

#' Evaluate SLD Composite Log-Density on New Data
#'
#' Applies the trained spectral shift to out-of-sample observations.
#'
#' @param newdata_encoded List with Z and Z_mis from
#'   \code{encode_newdata_spectral}.
#' @param encoding Training encoding list.
#' @param pi_c Flattened class-conditional marginal probabilities.
#' @param A_star Hollow projection matrix.
#' @return Numeric vector of length N (log-densities).
#' @keywords internal
eval_spectral_density_oos <- function(newdata_encoded, encoding, pi_c, A_star) {
  N <- nrow(newdata_encoded$Z)
  C <- encoding$C

  pi_c <- pmax(pi_c, 1e-15)
  E_c <- matrix(pi_c, nrow = N, ncol = C, byrow = TRUE)

  Z_imp <- newdata_encoded$Z
  Z_imp[newdata_encoded$Z_mis] <- E_c[newdata_encoded$Z_mis]

  R_c <- Z_imp - E_c
  S_c <- R_c %*% A_star
  eta <- log(E_c) + S_c
  eta <- pmin(pmax(eta, -30), 30)  # clamp for numerical stability
  ld <- numeric(N)

  for (v in names(encoding$item_indices)) {
    idx <- encoding$item_indices[[v]]
    eta_j <- eta[, idx, drop = FALSE]

    row_max <- apply(eta_j, 1, max, na.rm = TRUE)
    shifted <- exp(sweep(eta_j, 1, row_max, "-"))
    p_hat_j <- sweep(shifted, 1, rowSums(shifted, na.rm = TRUE), "/")
    p_hat_j <- pmax(p_hat_j, 1e-15)

    Z_j <- newdata_encoded$Z[, idx, drop = FALSE]
    log_P <- log(p_hat_j)
    ld <- ld + rowSums(Z_j * log_P)
  }
  ld
}

# ------------------------------------------------------------------------------
# Diagnostics and Visualization
# ------------------------------------------------------------------------------

#' Spectral Loadings as a Tidy Data Frame
#' @param model A \code{mixLCA} object.
#' @return Data frame.
#' @export
spectral_loadings <- function(model) {
  if (is.null(model$cat_spectral_params)) stop("`model` does not contain spectral local dependence parameters.")

  sp <- model$cat_spectral_params
  enc <- sp$encoding
  d_vec <- model$specs$spectral_rank
  if (length(d_vec) == 1L) d_vec <- rep(d_vec, model$n_classes)
  C <- enc$C

  col_item <- enc$item_names[enc$item_of]
  col_cat <- character(C)
  for (v in enc$item_names) col_cat[enc$item_indices[[v]]] <- enc$cat_levels[[v]]

  build_block <- function(V, class_label, d_k) {
    if (d_k < 1L) return(NULL)
    out <- expand.grid(col = seq_len(C), dimension = seq_len(d_k), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    out$class <- class_label
    out$item <- col_item[out$col]
    out$category <- col_cat[out$col]
    out$loading <- V[cbind(out$col, out$dimension)]
    out[, c("class", "dimension", "item", "category", "loading")]
  }

  if (isTRUE(model$specs$spectral_pool)) {
    df <- build_block(sp$V_d[[1L]], "Pooled", max(d_vec))
  } else {
    df <- do.call(rbind, lapply(seq_along(sp$V_d), function(k) build_block(sp$V_d[[k]], paste("Class", k), d_vec[k])))
  }
  df
}

#' Eigenvalue Scree Plot for a Spectral Model
#' @param model A \code{mixLCA} object.
#' @param max_eigs Maximum number of leading eigenvalues to display per class.
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_spectral_scree <- function(model, max_eigs = 20L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  if (is.null(model$cat_spectral_params)) stop("Model does not have spectral local dependence.")

  sp <- model$cat_spectral_params
  d_vec <- model$specs$spectral_rank
  if (length(d_vec) == 1L) d_vec <- rep(d_vec, model$n_classes)
  trim_eigs <- function(eigs, max_n) utils::head(eigs[abs(eigs) > 1e-10], max_n)

  if (isTRUE(model$specs$spectral_pool)) {
    eigs <- trim_eigs(sp$spectra[[1L]], max_eigs)
    df <- data.frame(index = seq_along(eigs), eigenvalue = eigs, class = "Pooled")
  } else {
    df <- do.call(rbind, lapply(seq_along(sp$spectra), function(k) {
      if (d_vec[k] == 0L) return(NULL)
      eigs <- trim_eigs(sp$spectra[[k]], max_eigs)
      data.frame(index = seq_along(eigs), eigenvalue = eigs, class = paste("Class", k))
    }))
  }

  d_line <- max(d_vec)
  ggplot2::ggplot(df, ggplot2::aes(x = index, y = eigenvalue, color = class, group = class)) +
    ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_vline(xintercept = d_line + 0.5, linetype = "dashed", color = "firebrick") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(x = "Eigenvalue index", y = "Eigenvalue", color = NULL,
                  title = "Conditional Burt Matrix Spectrum - mixLCA",
                  subtitle = paste0("Ranks d = [", paste(d_vec, collapse = ", "), "] (dashed at max)"))
}

#' Loadings Plot for a Spectral Dimension
#' @param model A \code{mixLCA} object.
#' @param dimension Integer in \code{1:d}.
#' @param class Integer class index.
#' @param n_top Optional integer: keep only top absolute loadings.
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_spectral_loadings <- function(model, dimension = 1L, class = NULL, n_top = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
  
  ld <- spectral_loadings(model)
  ld <- ld[ld$dimension == dimension, , drop = FALSE]

  if (isTRUE(model$specs$spectral_pool)) {
    df <- ld; title_suffix <- "Pooled"
  } else {
    if (is.null(class)) class <- 1L
    df <- ld[ld$class == paste("Class", class), , drop = FALSE]
    title_suffix <- paste("Class", class)
  }

  if (!is.null(n_top) && nrow(df) > n_top) df <- df[utils::head(order(abs(df$loading), decreasing = TRUE), n_top), ]

  df$column_label <- factor(paste(df$item, df$category, sep = " :: "), levels = paste(df$item, df$category, sep = " :: ")[order(df$loading)])

  ggplot2::ggplot(df, ggplot2::aes(x = loading, y = column_label, fill = item)) +
    ggplot2::geom_col(color = "grey30", linewidth = 0.3) + ggplot2::geom_vline(xintercept = 0, color = "grey30") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(x = paste("Loading on dimension", dimension), y = NULL, fill = "Item",
                  title = paste("Spectral Loadings:", title_suffix))
}
