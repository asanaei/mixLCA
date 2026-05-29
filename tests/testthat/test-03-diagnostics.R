# ==============================================================================
# tests/testthat/test-03-diagnostics.R
# Tests for fit_indices, bvr, bvr_tests, class_table,
# compare_models.
# ==============================================================================

test_that("fit_indices() returns all required statistics", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  fi <- fit_indices(fit)

  expected_names <- c("log_lik","n_params","AIC","BIC","aBIC",
                       "entropy","ICL","is_composite")
  expect_true(all(expected_names %in% names(fi)))
  expect_true(all(sapply(fi[expected_names], is.finite)))
})

test_that("fit_indices() entropy is in [0, 1]", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  fi <- fit_indices(fit)
  expect_gte(fi$entropy, 0)
  expect_lte(fi$entropy, 1)
})

test_that("fit_indices() BIC > AIC (more parameters penalized)", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  fi <- fit_indices(fit)
  # BIC penalty = np * log(N) >> AIC penalty = 2*np for N >> e^2
  expect_true(fi$BIC > fi$AIC)
})

test_that("bvr() returns a symmetric d x d matrix", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  bvr <- bvr(fit, df)

  expect_equal(nrow(bvr), 3L)
  expect_equal(ncol(bvr), 3L)
  expect_equal(rownames(bvr), c("x1","x2","x3"))
  expect_equal(max(abs(bvr - t(bvr))), 0)
})

test_that("bvr_tests() returns correct number of pairs", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  tests <- bvr_tests(fit, df)

  # 3 indicators -> 3 pairs
  expect_equal(nrow(tests), 3L)
  expect_true(all(c("var1","var2","residual_cov","chi_sq","p_value") %in%
                    names(tests)))
  expect_true(all(tests$p_value >= 0 & tests$p_value <= 1))
  expect_true(all(tests$chi_sq >= 0))
})

test_that("bvr() residuals are near zero when model fits well", {
  # Model with diagonal covariance on data generated under independence
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 dependence = "none",
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  bvr <- bvr(fit, df)

  # Off-diagonal residuals should be modest (data have zero true off-diag cov)
  off <- bvr[upper.tri(bvr)]
  expect_true(all(abs(off) < 5),
              label = "off-diagonal residuals are not extreme")
})

test_that("class_table() returns K x K matrix with row sums near 1", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  ct <- class_table(fit)

  expect_equal(dim(ct), c(2L, 2L))
  expect_true(all(abs(rowSums(ct) - 1) < 1e-6))
  expect_true(all(ct >= 0))
})

test_that("compare_models() returns a data frame ordered by K", {
  df  <- make_df_small()
  f2  <- fit_lca(df, continuous = c("x1","x2","x3"), n_classes = 2,
                 control = lca_control(n_starts = 1), verbose = FALSE)
  f3  <- fit_lca(df, continuous = c("x1","x2","x3"), n_classes = 3,
                 control = lca_control(n_starts = 1), verbose = FALSE)
  tab <- compare_models(f2, f3)

  expect_s3_class(tab, "data.frame")
  expect_equal(nrow(tab), 2L)
  expect_true(all(diff(tab$K) >= 0))
  expect_true("BIC" %in% names(tab))
})

test_that("fit_indices() errors on non-mixLCA input", {
  expect_error(fit_indices(list()), regexp = "mixLCA")
})
