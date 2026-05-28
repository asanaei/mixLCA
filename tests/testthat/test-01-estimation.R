# ==============================================================================
# tests/testthat/test-01-estimation.R
# Tests for fit_lca() — class structure, field types, parameter recovery.
# ==============================================================================

test_that("fit_lca() returns a mixLCA object with required fields", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)

  expect_s3_class(fit, "mixLCA")

  required <- c("n_classes", "log_lik", "posteriors",
                 "continuous_params", "n_obs", "n_params",
                 "convergence", "specs")
  for (fld in required)
    expect_true(!is.null(fit[[fld]]),
                label = paste("field present:", fld))
})

test_that("fit_lca() posterior matrix has correct shape and row sums", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)

  post <- fit$posteriors
  expect_equal(nrow(post), nrow(df))
  expect_equal(ncol(post), 2L)
  expect_true(all(post >= 0))
  expect_true(all(post <= 1 + 1e-9))
  expect_true(all(abs(rowSums(post) - 1) < 1e-8))
})

test_that("fit_lca() log-likelihood is a finite scalar", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)

  expect_true(is.finite(fit$log_lik))
  expect_true(fit$log_lik < 0)
})

test_that("fit_lca() n_params is a positive integer", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)

  expect_true(is.integer(fit$n_params))
  expect_true(fit$n_params > 0L)
})

test_that("fit_lca() ll_history is non-decreasing", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)

  ll <- fit$convergence$ll_history
  diffs <- diff(ll)
  # Allow tiny numerical noise (< 1e-4) but no substantial decreases
  expect_true(all(diffs > -1e-4),
              label = "log-likelihood is non-decreasing across EM iterations")
})

test_that("fit_lca() works with only categorical indicators", {
  df <- make_df_small()
  expect_no_error(
    fit_lca(df,
            categorical = "cat1",
            n_classes   = 2,
            control = lca_control(n_starts = 1),
            verbose     = FALSE)
  )
})

test_that("fit_lca() works with mixed continuous + categorical indicators", {
  df <- make_df_small()
  expect_no_error(
    fit_lca(df,
            continuous  = c("x1","x2"),
            categorical = "cat1",
            n_classes   = 2,
            control = lca_control(n_starts = 1),
            verbose     = FALSE)
  )
})

test_that("fit_lca() works with concomitant predictors", {
  df <- make_df_small()
  fit <- fit_lca(df,
                 continuous  = c("x1","x2","x3"),
                 concomitant = "age",
                 n_classes   = 2,
                 control = lca_control(n_starts = 1),
                 verbose     = FALSE)

  expect_false(is.null(fit$concomitant_coefs))
  expect_equal(dim(fit$concomitant_coefs), c(2L, 1L))  # P=2 (intercept + age), K-1=1
})

test_that("fit_lca() handles missing continuous data without error", {
  df      <- make_df_small()
  df$x1[sample(nrow(df), 30)] <- NA

  expect_no_error(
    fit_lca(df,
            continuous = c("x1","x2","x3"),
            n_classes  = 2,
            control = lca_control(n_starts = 1),
            verbose    = FALSE)
  )
})

test_that("fit_lca() with dependence='none' produces diagonal covariances", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 dependence = "none",
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)

  for (k in seq_len(fit$n_classes)) {
    Sk <- fit$continuous_params$covariances[[k]]
    off_diag <- Sk[upper.tri(Sk)]
    expect_true(all(abs(off_diag) < 1e-9),
                label = paste("class", k, "covariance is diagonal"))
  }
})

test_that("fit_lca() with dependence='penalized' returns non-negative variances", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 dependence = "penalized",
                 penalty    = 1.0,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)

  for (k in seq_len(fit$n_classes))
    expect_true(all(diag(fit$continuous_params$covariances[[k]]) > 0))
})

test_that("fit_lca() class means differ between classes for well-separated data", {
  df  <- make_df_small()   # x1: class1 ~ N(10,2), class2 ~ N(4,2)
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 2),
                 verbose    = FALSE)

  m1 <- fit$continuous_params$means[[1]]
  m2 <- fit$continuous_params$means[[2]]
  # At least one indicator should have means that differ by > 2
  max_sep <- max(abs(m1 - m2))
  expect_true(max_sep > 2,
              label = "classes are meaningfully separated")
})

test_that("multi-start produces equal or higher LL than single start", {
  df   <- make_df_small()
  fit1 <- fit_lca(df,
                  continuous = c("x1","x2","x3"),
                  n_classes  = 2,
                  control = lca_control(n_starts = 1),
                  verbose    = FALSE)
  fit5 <- fit_lca(df,
                  continuous = c("x1","x2","x3"),
                  n_classes  = 2,
                  control = lca_control(n_starts = 5),
                  verbose    = FALSE)

  expect_true(fit5$log_lik >= fit1$log_lik - 1e-4)
})

test_that("validate_inputs catches bad inputs", {
  df <- make_df_small()
  expect_error(fit_lca(df, continuous = c("x1","x2","x3"),
                        n_classes = 0, verbose = FALSE),
               regexp = "n_classes")
  expect_error(fit_lca(df, continuous = "nonexistent",
                        n_classes = 2, verbose = FALSE),
               regexp = "absent")
  expect_error(fit_lca(df, n_classes = 2, verbose = FALSE),
               regexp = "least one")
})
