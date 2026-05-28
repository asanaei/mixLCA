# ==============================================================================
# tests/testthat/test-04-distal.R
# Tests for distal() and the mixDistal S3 class.
# ==============================================================================

test_that("distal() returns a mixDistal object with required fields", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  dis <- distal(fit, df, outcome ~ 1, "gaussian")

  expect_s3_class(dis, "mixDistal")
  expect_true(!is.null(dis$class_models))
  expect_equal(length(dis$class_models), 2L)
  expect_equal(dis$n_classes, 2L)
})

test_that("distal() unconditional estimates differ between classes", {
  df  <- make_df_small()  # outcome: class1 ~ N(80,10), class2 ~ N(40,10)
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 2),
                 verbose    = FALSE)
  dis <- distal(fit, df, outcome ~ 1, "gaussian")

  mu1 <- dis$class_models[[1]]$coef_table["(Intercept)", "Estimate"]
  mu2 <- dis$class_models[[2]]$coef_table["(Intercept)", "Estimate"]
  # The two means should be at least 15 units apart
  expect_true(abs(mu1 - mu2) > 15,
              label = "class-specific distal means are well separated")
})

test_that("distal() with predictor produces non-trivial slope", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 2),
                 verbose    = FALSE)
  dis <- distal(fit, df, outcome ~ age, "gaussian")

  for (k in 1:2) {
    tab <- dis$class_models[[k]]$coef_table
    expect_true(all(is.finite(tab)),
                label = paste("all estimates finite for class", k))
  }
})

test_that("distal() classification error matrix columns sum to 1", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  dis <- distal(fit, df, outcome ~ 1, "gaussian")

  col_sums <- colSums(dis$error_matrix)
  expect_true(all(abs(col_sums - 1) < 1e-6))
})

test_that("distal() errors on non-mixLCA input", {
  expect_error(distal(list(), data.frame(), ~1),
               regexp = "mixLCA")
})

test_that("print.mixDistal produces output without error", {
  df  <- make_df_small()
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  dis <- distal(fit, df, outcome ~ 1, "gaussian")

  expect_output(print(dis), "BCH")
  expect_output(print(dis), "Class 1")
  expect_output(print(dis), "Class 2")
})

test_that("distal() handles missing outcome values", {
  df               <- make_df_small()
  df$outcome[1:10] <- NA
  fit <- fit_lca(df,
                 continuous = c("x1","x2","x3"),
                 n_classes  = 2,
                 control = lca_control(n_starts = 1),
                 verbose    = FALSE)
  expect_no_error(distal(fit, df, outcome ~ 1, "gaussian"))
})
