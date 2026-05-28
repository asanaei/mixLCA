# ==============================================================================
# tests/testthat/test-02-s3-methods.R
# Tests for print, summary, predict, coef, logLik, nobs.
# ==============================================================================

# Shared fixture for this file
fit_small <- function() {
  fit_lca(make_df_small(),
          continuous  = c("x1","x2","x3"),
          concomitant = "age",
          n_classes   = 2,
          control = lca_control(n_starts = 1),
          verbose     = FALSE)
}

test_that("print.mixLCA produces output without error", {
  fit <- fit_small()
  expect_output(print(fit), "mixLCA")
  expect_output(print(fit), "Classes")
  expect_output(print(fit), "Log-likelihood")
})

test_that("summary.mixLCA produces output without error", {
  fit <- fit_small()
  expect_output(summary(fit), "mixLCA")
  expect_output(summary(fit), "AIC")
  expect_output(summary(fit), "BIC")
  expect_output(summary(fit), "Entropy")
})

test_that("summary.mixLCA includes SE table when data is supplied", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous  = c("x1","x2","x3"),
                 concomitant = "age",
                 n_classes   = 2,
                 control = lca_control(n_starts = 1),
                 verbose     = FALSE)
  expect_output(summary(fit, data = df), "Estimate")
})

test_that("predict.mixLCA default returns posterior probability matrix", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  prob <- predict(fit)

  expect_true(is.matrix(prob))
  expect_equal(nrow(prob), nrow(df))
  expect_equal(ncol(prob), 2L)
  expect_true(all(abs(rowSums(prob) - 1) < 1e-8))
})

test_that("predict.mixLCA type = 'class' returns modal integer vector", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  cls <- predict(fit, type = "class")

  expect_type(cls, "integer")
  expect_length(cls, nrow(df))
  expect_true(all(cls %in% 1:2))
})

test_that("predict.mixLCA type = 'all' returns the legacy data frame", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  preds <- predict(fit, type = "all")

  expect_s3_class(preds, "data.frame")
  expect_equal(nrow(preds), nrow(df))
  expect_equal(ncol(preds), 4L)
  expect_true(all(preds$max_posterior >= 0.5 - 1e-9))
  expect_true(all(preds$max_posterior <= 1 + 1e-9))
  expect_true(all(preds$modal_class %in% 1:2))
})

test_that("coef.mixLCA returns a named list with expected elements", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous  = c("x1","x2","x3"),
                 concomitant = "age",
                 n_classes   = 2,
                 control = lca_control(n_starts = 1),
                 verbose     = FALSE)
  cc <- coef(fit)

  expect_type(cc, "list")
  expect_true("concomitant" %in% names(cc))
  expect_true("means" %in% names(cc))
  expect_true("covariances" %in% names(cc))
})

test_that("logLik.mixLCA returns correct class and attributes", {
  fit <- fit_small()
  ll  <- logLik(fit)

  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "df"),   fit$n_params)
  expect_equal(attr(ll, "nobs"), fit$n_obs)
  expect_equal(as.numeric(ll),   fit$log_lik)
})

test_that("nobs.mixLCA returns the correct count", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  expect_equal(nobs(fit), nrow(df))
})
