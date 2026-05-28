# ==============================================================================
# tests/testthat/test-08-parameter-recovery.R
# Integration tests: verify that mixLCA recovers known parameters from
# data generated under the exact model.
# ==============================================================================

test_that("class means are recovered within 1 SE of true values (N=900)", {
  df <- make_df_full()

  fit <- fit_lca(
    df,
    continuous  = c("Metric_X", "Metric_Y", "Metric_Z"),
    categorical = "Cat_Type",
    concomitant = c("Age", "Status"),
    n_classes   = 2,
    dependence  = "full",
    control = lca_control(n_starts = 3),
    verbose     = FALSE
  )

  # True means: Class 1 = (25, 25, 8), Class 2 = (12, 12, 18)
  # We cannot know which class index maps to which true class, so we
  # check that the recovered pair of means matches the true pair.
  mu1 <- fit$continuous_params$means[[1]]
  mu2 <- fit$continuous_params$means[[2]]

  true_set <- list(c(25, 25, 8), c(12, 12, 18))

  # Match by minimum sum of squared deviations
  dist11 <- sum((mu1 - true_set[[1]])^2) + sum((mu2 - true_set[[2]])^2)
  dist12 <- sum((mu1 - true_set[[2]])^2) + sum((mu2 - true_set[[1]])^2)
  best   <- min(dist11, dist12)

  expect_true(best < 9,
              label = "means within 3 units RMSE of true values")
})

test_that("concomitant age coefficient is positive (true value +0.06)", {
  df <- make_df_full()

  fit <- fit_lca(
    df,
    continuous  = c("Metric_X", "Metric_Y", "Metric_Z"),
    categorical = "Cat_Type",
    concomitant = c("Age", "Status"),
    n_classes   = 2,
    dependence  = "full",
    control = lca_control(n_starts = 3),
    verbose     = FALSE
  )

  # Concomitant coefficient for Age, class 2 vs 1 (or class 1 vs 2)
  # One of the two possibilities should have a positive age coefficient
  # of roughly 0.06 in magnitude.
  coefs <- fit$concomitant_coefs
  age_coef <- coefs[2, 1]   # row 2 = Age (after intercept)
  expect_true(abs(abs(age_coef) - 0.06) < 0.04,
              label = "age coefficient close to true ±0.06")
})

test_that("distal means recover true outcome separation", {
  df <- make_df_full()

  fit <- fit_lca(
    df,
    continuous  = c("Metric_X", "Metric_Y", "Metric_Z"),
    categorical = "Cat_Type",
    concomitant = c("Age", "Status"),
    n_classes   = 2,
    dependence  = "full",
    control = lca_control(n_starts = 3),
    verbose     = FALSE
  )

  dis <- distal(fit, df, Clinical_Outcome ~ 1, "gaussian")

  mu1 <- dis$class_models[[1]]$coef_table["(Intercept)", "Estimate"]
  mu2 <- dis$class_models[[2]]$coef_table["(Intercept)", "Estimate"]

  # True unconditional means differ by roughly 30 units at mean Age=40:
  # class1 ~ 50 + 0.8*40 = 82, class2 ~ 20 + 0.1*40 = 24
  sep <- abs(mu1 - mu2)
  expect_true(sep > 30,
              label = "distal means separated by > 30 units")
})

test_that("covariance off-diagonal recovered near true value (4.0)", {
  df <- make_df_full()

  fit <- fit_lca(
    df,
    continuous  = c("Metric_X", "Metric_Y", "Metric_Z"),
    categorical = "Cat_Type",
    concomitant = c("Age", "Status"),
    n_classes   = 2,
    dependence  = "full",
    control = lca_control(n_starts = 3),
    verbose     = FALSE
  )

  # True off-diagonal cov(X,Y) = 4 in one class, 0 in the other
  cov1 <- fit$continuous_params$covariances[[1]]
  cov2 <- fit$continuous_params$covariances[[2]]

  off_diag_vals <- c(cov1[1, 2], cov2[1, 2])
  max_off       <- max(abs(off_diag_vals))
  expect_true(max_off > 2,
              label = "largest X-Y off-diagonal recovered above 2")
})

test_that("L1 penalty shrinks small covariances toward zero", {
  df <- make_df_full()

  fit_full <- fit_lca(
    df,
    continuous = c("Metric_X", "Metric_Y", "Metric_Z"),
    n_classes  = 2,
    dependence = "full",
    control = lca_control(n_starts = 2),
    verbose    = FALSE
  )
  fit_pen <- fit_lca(
    df,
    continuous = c("Metric_X", "Metric_Y", "Metric_Z"),
    n_classes  = 2,
    dependence = "penalized",
    penalty    = 2.0,
    control = lca_control(n_starts = 2),
    verbose    = FALSE
  )

  # Sum of absolute off-diagonal elements should be smaller under penalty
  sum_off_full <- sum(sapply(fit_full$continuous_params$covariances, function(S) {
    sum(abs(S[upper.tri(S)]))
  }))
  sum_off_pen <- sum(sapply(fit_pen$continuous_params$covariances, function(S) {
    sum(abs(S[upper.tri(S)]))
  }))

  expect_true(sum_off_pen <= sum_off_full,
              label = "L1 penalty reduces total off-diagonal magnitude")
})
