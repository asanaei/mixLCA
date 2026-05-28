# ==============================================================================
# File: R/09_visualization.R
# mixLCA - ggplot2-based visualisations.
# ==============================================================================

#' Plot Method for mixLCA
#'
#' Single dispatch surface for all built-in visualisations. The
#' \code{type} argument selects which plot to produce.
#'
#' Plot-specific arguments (e.g.\ \code{data} for \code{"bvr"} and
#' \code{"distal"}, \code{variable} for \code{"distal"}, \code{ci}
#' for \code{"profiles"}, \code{dimension} for
#' \code{"spectral_loadings"}) are forwarded via \code{...}.
#'
#' @param x A fitted \code{mixLCA} object.
#' @param type One of \code{"profiles"}, \code{"bvr"}, \code{"distal"},
#'   \code{"uncertainty"}, \code{"convergence"}, \code{"categorical"},
#'   \code{"spectral_scree"}, or \code{"spectral_loadings"}.
#' @param ... Forwarded to the underlying plotting routine.
#' @return A \code{ggplot} object.
#' @export
plot.mixLCA <- function(x, type = c("profiles", "bvr", "distal",
                                    "uncertainty", "convergence",
                                    "categorical", "spectral_scree",
                                    "spectral_loadings"), ...) {
  type <- match.arg(type)
  switch(type,
    profiles          = .plot_profiles(x, ...),
    bvr               = .plot_bvr(x, ...),
    distal            = .plot_distal(x, ...),
    uncertainty       = .plot_uncertainty(x, ...),
    convergence       = .plot_convergence(x, ...),
    categorical       = .plot_categorical(x, ...),
    spectral_scree    = .plot_spectral_scree(x, ...),
    spectral_loadings = .plot_spectral_loadings(x, ...)
  )
}

#' Profile Plot of Continuous Indicator Means
#'
#' Renders class-specific mean profiles with optional 95\% confidence
#' intervals derived from numerical standard errors. Classes are
#' optionally reordered by profile severity (sum of means) to produce
#' a consistent severity gradient across model runs.
#'
#' @param model A \code{mixLCA} object with continuous indicators.
#' @param data Data frame (required if \code{ci = TRUE}).
#' @param ci Logical: draw 95\% confidence intervals?
#' @param reorder Logical: reorder classes by profile severity?
#' @param ... Additional arguments passed to
#'   \code{\link[ggplot2]{geom_line}} (e.g. \code{linewidth}).
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_profiles <- function(model, data = NULL, ci = FALSE,
                           reorder = TRUE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plot(type = 'profiles').")
  if (!inherits(model, "mixLCA"))
    stop("`model` must be a mixLCA object.")

  cont  <- model$specs$continuous
  if (is.null(cont)) stop("No continuous indicators to plot.")

  K     <- model$n_classes
  props <- colMeans(model$posteriors)
  mu_list <- model$continuous_params$means

  ord <- if (reorder) order(sapply(mu_list, sum)) else seq_len(K)

  se_list <- NULL
  if (ci && !is.null(data)) {
    cse     <- continuous_se(model, data)
    se_list <- cse$mean_se
  }

  frames <- list()
  for (rank in seq_len(K)) {
    k     <- ord[rank]
    label <- paste0("Class ", rank, " (",
                    round(100 * props[k], 1), "%)")
    fr <- data.frame(
      indicator = factor(cont, levels = cont),
      mean      = mu_list[[k]],
      class     = label,
      stringsAsFactors = FALSE
    )
    if (!is.null(se_list)) {
      fr$lower <- mu_list[[k]] - 1.96 * se_list[[k]]
      fr$upper <- mu_list[[k]] + 1.96 * se_list[[k]]
    }
    frames[[rank]] <- fr
  }
  plot_df <- do.call(rbind, frames)

  p <- ggplot2::ggplot(plot_df,
         ggplot2::aes(x = indicator, y = mean,
                      group = class, color = class)) +
    ggplot2::geom_line(linewidth = 1, ...) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(x = NULL, y = "Mean", color = NULL,
                  title = "Measurement Profiles - mixLCA")

  if (!is.null(se_list))
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = lower, ymax = upper),
        width = 0.15, linewidth = 0.5)
  p
}

#' Bivariate Residual Network Graph
#'
#' Constructs a network graph where nodes are continuous indicators
#' and edges connect pairs whose residual covariance is significant
#' at p < .05 (default) or exceeds a user-supplied numeric threshold.
#' Edge thickness and opacity encode residual magnitude; edge colour
#' distinguishes positive (red) from negative (blue) residuals.
#'
#' @param model A \code{mixLCA} object.
#' @param data Data frame used for estimation.
#' @param threshold Numeric minimum absolute residual covariance, or
#'   \code{"sig"} to filter by chi-squared p < .05.
#' @return A \code{ggplot} object via \code{ggraph}, or NULL when no
#'   pairs exceed the threshold.
#' @keywords internal
.plot_bvr <- function(model, data, threshold = "sig") {
  if (!requireNamespace("igraph",  quietly = TRUE))
    stop("igraph is required for plot(type = 'bvr').")
  if (!requireNamespace("ggraph",  quietly = TRUE))
    stop("ggraph is required for plot(type = 'bvr').")
  if (!inherits(model, "mixLCA"))
    stop("`model` must be a mixLCA object.")

  tests <- bvr_tests(model, data)

  edges <- if (is.character(threshold) && threshold == "sig") {
    tests[tests$p_value < 0.05, ]
  } else {
    tests[abs(tests$residual_cov) > threshold, ]
  }

  if (nrow(edges) == 0L) {
    message("No bivariate residuals exceed the threshold. ",
            "Local independence holds adequately.")
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = "Bivariate Residual Network - mixLCA",
                           subtitle = "No edges exceed the threshold (local independence holds)."))
  }

  g <- igraph::graph_from_data_frame(
    edges[, c("var1", "var2", "residual_cov")],
    directed = FALSE
  )
  igraph::E(g)$magnitude <- abs(edges$residual_cov)
  igraph::E(g)$sign      <- ifelse(edges$residual_cov > 0,
                                   "positive", "negative")

  ggraph::ggraph(g, layout = "circle") +
    ggraph::geom_edge_link(
      ggplot2::aes(width = magnitude, alpha = magnitude,
                   color = sign)) +
    ggraph::scale_edge_color_manual(
      values = c(positive = "firebrick", negative = "steelblue")) +
    ggraph::scale_edge_width(range = c(0.5, 3)) +
    ggraph::geom_node_point(size = 8, color = "grey30") +
    ggraph::geom_node_text(ggplot2::aes(label = name),
                           repel = TRUE, size = 4.5) +
    ggplot2::theme_void(base_size = 13) +
    ggplot2::labs(
      title    = "Bivariate Residual Network - mixLCA",
      subtitle = "Edges denote local dependence violations")
}

#' Distal Outcome Density by Class
#'
#' Plots kernel density estimates of a distal variable split by modal
#' class assignment. Observations may be weighted by their maximum
#' posterior probability to reflect classification uncertainty.
#'
#' @param model A \code{mixLCA} object.
#' @param data Data frame containing the distal variable.
#' @param variable Character: name of the distal variable.
#' @param weighted Logical: weight densities by max posterior
#'   probability?
#' @param ... Additional arguments passed to
#'   \code{\link[ggplot2]{geom_density}} (e.g. \code{bw}, \code{adjust},
#'   \code{kernel}).
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_distal <- function(model, data, variable, weighted = TRUE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plot(type = 'distal').")
  if (!inherits(model, "mixLCA"))
    stop("`model` must be a mixLCA object.")

  modal  <- apply(model$posteriors, 1, which.max)
  max_p  <- apply(model$posteriors, 1, max)
  plot_df <- data.frame(
    value  = data[[variable]],
    class  = paste("Class", modal),
    weight = if (weighted) max_p else 1,
    stringsAsFactors = FALSE
  )
  plot_df <- plot_df[!is.na(plot_df$value), ]

  p <- if (weighted) {
    ggplot2::ggplot(plot_df,
      ggplot2::aes(x = value, weight = weight, fill = class)) +
    ggplot2::geom_density(alpha = 0.45, color = "grey30", ...) +
    ggplot2::labs(subtitle = "Weighted by classification certainty")
  } else {
    ggplot2::ggplot(plot_df,
      ggplot2::aes(x = value, fill = class)) +
    ggplot2::geom_density(alpha = 0.45, color = "grey30", ...)
  }

  p + ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(x = variable, y = "Density", fill = NULL,
                  title = paste("Distal Distribution:", variable,
                                "- mixLCA"))
}

#' Classification Uncertainty Histogram
#'
#' Displays the distribution of maximum posterior probabilities across
#' observations. A peak near 1.0 indicates confident classification;
#' mass near 1/K indicates substantial ambiguity.
#'
#' @param model A \code{mixLCA} object.
#' @param ... Additional arguments passed to
#'   \code{\link[ggplot2]{geom_histogram}} (e.g. \code{bins},
#'   \code{binwidth}).
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_uncertainty <- function(model, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plot(type = 'uncertainty').")
  if (!inherits(model, "mixLCA"))
    stop("`model` must be a mixLCA object.")

  max_p   <- apply(model$posteriors, 1, max)
  plot_df <- data.frame(max_posterior = max_p)

  ggplot2::ggplot(plot_df, ggplot2::aes(x = max_posterior)) +
    ggplot2::geom_histogram(fill = "grey50",
                            color = "white", linewidth = 0.3, ...) +
    ggplot2::geom_vline(xintercept = 1 / model$n_classes,
                        linetype = "dashed", color = "firebrick") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(
      x        = "Maximum posterior probability",
      y        = "Count",
      title    = "Classification Certainty - mixLCA",
      subtitle = paste("Dashed line = chance level (1 /",
                       model$n_classes, ")"))
}

#' EM Convergence Trace Plot
#'
#' Plots the log-likelihood across EM iterations.
#'
#' @param model A \code{mixLCA} object.
#' @param ... Additional arguments passed to
#'   \code{\link[ggplot2]{geom_line}}.
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_convergence <- function(model, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plot(type = 'convergence').")
  if (!inherits(model, "mixLCA"))
    stop("`model` must be a mixLCA object.")

  ll      <- model$convergence$ll_history
  plot_df <- data.frame(iteration = seq_along(ll), log_lik = ll)

  ggplot2::ggplot(plot_df,
    ggplot2::aes(x = iteration, y = log_lik)) +
    ggplot2::geom_line(linewidth = 0.8, ...) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(x = "Iteration", y = "Log-likelihood",
                  title = "EM Convergence - mixLCA")
}

#' Plot Categorical Indicator Probabilities
#'
#' Visualizes the class-conditional response probabilities for categorical
#' manifest variables using dodged bar charts faceted by variable.
#'
#' @param model A \code{mixLCA} object.
#' @param variables Character vector: optional subset of categorical
#'   variables to plot. If \code{NULL}, plots all categorical indicators.
#' @param orientation Character: \code{"vertical"} (default) or
#'   \code{"horizontal"}. Horizontal is recommended for variables with
#'   long category names.
#' @return A \code{ggplot} object.
#' @keywords internal
.plot_categorical <- function(model, variables = NULL,
                              orientation = c("vertical", "horizontal")) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plot(type = 'categorical').")
  if (!inherits(model, "mixLCA"))
    stop("`model` must be a mixLCA object.")

  cat_params <- model$categorical_params
  if (is.null(cat_params))
    stop("No categorical parameters found in the model.")

  orientation <- match.arg(orientation)
  K <- model$n_classes

  all_vars <- names(cat_params[[1]])
  if (is.null(variables)) {
    variables <- all_vars
  } else {
    missing_vars <- setdiff(variables, all_vars)
    if (length(missing_vars) > 0L)
      stop("Variables not found in categorical parameters: ",
           paste(missing_vars, collapse = ", "))
  }

  # Flatten nested parameter list into a tidy plotting dataframe
  plot_rows <- list()
  idx <- 1L
  for (k in seq_len(K)) {
    for (v in variables) {
      probs <- cat_params[[k]][[v]]
      for (cat_name in names(probs)) {
        plot_rows[[idx]] <- data.frame(
          Class       = paste("Class", k),
          Variable    = v,
          Category    = as.character(cat_name),
          Probability = probs[[cat_name]],
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  df <- do.call(rbind, plot_rows)

  df$Class    <- factor(df$Class,    levels = paste("Class", seq_len(K)))
  df$Variable <- factor(df$Variable, levels = variables)
  df$Category <- factor(df$Category, levels = unique(df$Category))

  p <- ggplot2::ggplot(df)
  if (orientation == "horizontal") {
    p <- p +
      ggplot2::aes(x = Probability, y = Category, fill = Class) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                        width = 0.7, color = "grey30", linewidth = 0.3) +
      ggplot2::facet_wrap(~ Variable, scales = "free_y") +
      ggplot2::scale_x_continuous(limits = c(0, 1),
                                  labels = function(x) paste0(x * 100, "%")) +
      ggplot2::labs(x = "Response Probability", y = "Category")
  } else {
    p <- p +
      ggplot2::aes(x = Category, y = Probability, fill = Class) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                        width = 0.7, color = "grey30", linewidth = 0.3) +
      ggplot2::facet_wrap(~ Variable, scales = "free_x") +
      ggplot2::scale_y_continuous(limits = c(0, 1),
                                  labels = function(x) paste0(x * 100, "%")) +
      ggplot2::labs(x = "Category", y = "Response Probability")
  }

  p + ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(fill = "Latent Class",
                  title = "Categorical Item Probabilities - mixLCA") +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"))
}
