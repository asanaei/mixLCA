# ==============================================================================
# tests/testthat/test-05-standard-errors.R
# Tests for concomitant_se() and continuous_se().
# ==============================================================================

test_that("concomitant_se() returns a matrix of the correct shape", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous  = c("x1","x2","x3"),
                 concomitant = "age",
                 n_classes   = 2,
                 control = lca_control(n_starts = 1),
                 verbose     = FALSE)
  se <- concomitant_se(fit, df)

  # P = 2 (intercept + age), K-1 = 1
  expect_equal(dim(se), c(2L, 1L))
  expect_true(all(is.finite(se)))
  expect_true(all(se > 0))
})

test_that("concomitant_se() returns NULL without concomitant predictors", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  expect_null(concomitant_se(fit, df))
})

test_that("concomitant_se() SEs are plausible relative to estimates", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous  = c("x1","x2","x3"),
                 concomitant = "age",
                 n_classes   = 2,
                 control = lca_control(n_starts = 2),
                 verbose     = FALSE)
  se <- concomitant_se(fit, df)
  ce <- fit$concomitant_coefs

  # z-statistics should be finite
  z <- ce / se
  expect_true(all(is.finite(z)))
})

test_that("continuous_se() returns lists with correct lengths", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  cse <- continuous_se(fit, df)

  expect_equal(length(cse$mean_se), 2L)
  expect_equal(length(cse$cov_se),  2L)
  expect_equal(length(cse$mean_se[[1]]), 3L)
})

test_that("continuous_se() SEs are positive or NA (not negative)", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  cse <- continuous_se(fit, df)

  for (k in 1:2) {
    se_k <- cse$mean_se[[k]]
    non_na <- se_k[!is.na(se_k)]
    expect_true(all(non_na > 0),
                label = paste("class", k, "mean SEs positive"))
  }
})

test_that("continuous_se() returns NULL without continuous indicators", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 categorical = "cat1",
                 n_classes   = 2,
                 control = lca_control(n_starts = 1),
                 verbose     = FALSE)
  expect_null(continuous_se(fit, df))
})

test_that("concomitant_se() errors on non-mixLCA input", {
  expect_error(concomitant_se(list(), data.frame()),
               regexp = "mixLCA")
})
