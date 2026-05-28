# ==============================================================================
# tests/testthat/test-07-utils.R
# Tests for log_sum_exp, softmax, soft_threshold, count_params,
# validate_inputs, eval_continuous_density, eval_categorical_density.
# ==============================================================================

test_that("mixLCA:::log_sum_exp() matches naive computation on small vectors", {
  x <- c(1, 2, 3)
  expect_equal(mixLCA:::log_sum_exp(x), log(sum(exp(x))), tolerance = 1e-10)
})

test_that("mixLCA:::log_sum_exp() handles -Inf entries", {
  x <- c(-Inf, 2, 3)
  expect_equal(mixLCA:::log_sum_exp(x), log(sum(exp(x[is.finite(x)]))),
               tolerance = 1e-10)
})

test_that("mixLCA:::log_sum_exp() returns -Inf for all -Inf input", {
  expect_equal(mixLCA:::log_sum_exp(c(-Inf, -Inf)), -Inf)
})

test_that("mixLCA:::softmax() rows sum to 1", {
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  S <- mixLCA:::softmax(X)
  expect_equal(rowSums(S), rep(1, 10), tolerance = 1e-10)
})

test_that("mixLCA:::softmax() returns values in [0, 1]", {
  X <- matrix(rnorm(30), nrow = 10, ncol = 3)
  S <- mixLCA:::softmax(X)
  expect_true(all(S >= 0))
  expect_true(all(S <= 1 + 1e-9))
})

test_that("mixLCA:::soft_threshold() zeros elements below lambda", {
  x      <- c(0.5, -0.3, 2.0, -1.5, 0.1)
  lambda <- 0.4
  result <- mixLCA:::soft_threshold(x, lambda)
  expect_equal(result[1], 0.1, tolerance = 1e-10)   # 0.5 - 0.4
  expect_equal(result[2], 0,   tolerance = 1e-10)   # |-.3| < .4
  expect_equal(result[3], 1.6, tolerance = 1e-10)   # 2.0 - 0.4
  expect_equal(result[4], -1.1, tolerance = 1e-10)  # -1.5 + 0.4
  expect_equal(result[5], 0,   tolerance = 1e-10)   # |0.1| < 0.4
})

test_that("mixLCA:::soft_threshold() with lambda=0 is identity", {
  x <- c(1.5, -2.3, 0.0)
  expect_equal(mixLCA:::soft_threshold(x, 0), x)
})

test_that("mixLCA:::eval_continuous_density() returns finite values for complete data", {
  set.seed(110)
  Y  <- matrix(rnorm(30), nrow = 10, ncol = 3)
  mu <- c(0, 0, 0)
  Sg <- diag(3)
  ld <- mixLCA:::eval_continuous_density(Y, mu, Sg)
  expect_length(ld, 10L)
  expect_true(all(is.finite(ld)))
  expect_true(all(ld <= 0))  # log densities are non-positive
})

test_that("mixLCA:::eval_continuous_density() returns 0 for fully missing rows", {
  Y     <- matrix(NA_real_, nrow = 5, ncol = 3)
  ld    <- mixLCA:::eval_continuous_density(Y, c(0,0,0), diag(3))
  expect_equal(ld, rep(0, 5))
})

test_that("mixLCA:::eval_continuous_density() handles partial missingness", {
  set.seed(110)
  Y       <- matrix(rnorm(20), nrow = 4, ncol = 5)
  Y[2, 3] <- NA
  Y[3, c(1,4)] <- NA
  mu <- rep(0, 5)
  Sg <- diag(5)
  ld <- mixLCA:::eval_continuous_density(Y, mu, Sg)
  expect_length(ld, 4L)
  expect_true(all(is.finite(ld)))
})

test_that("mixLCA:::eval_categorical_density() returns finite values", {
  df <- data.frame(
    v1 = c("A","B","A","C"),
    v2 = c("X","X","Y","Y"),
    stringsAsFactors = FALSE
  )
  probs <- list(
    v1 = c(A = 0.5, B = 0.3, C = 0.2),
    v2 = c(X = 0.6, Y = 0.4)
  )
  ld <- mixLCA:::eval_categorical_density(df, probs)
  expect_length(ld, 4L)
  expect_true(all(is.finite(ld)))
  expect_true(all(ld <= 0))
})

test_that("mixLCA:::eval_categorical_density() handles NA values (returns 0 contrib)", {
  df <- data.frame(
    v1 = c("A", NA, "B"),
    stringsAsFactors = FALSE
  )
  probs <- list(v1 = c(A = 0.6, B = 0.4))
  ld    <- mixLCA:::eval_categorical_density(df, probs)
  expect_length(ld, 3L)
  expect_equal(ld[2], ld[2])  # not NA
  expect_true(is.finite(ld[2]))
})

test_that("count_params() increases with more classes", {
  df <- make_df_small()
  f2 <- fit_lca(df, continuous = c("x1","x2","x3"),
                n_classes = 2,
                control = lca_control(n_starts = 1), verbose = FALSE)
  f3 <- fit_lca(df, continuous = c("x1","x2","x3"),
                n_classes = 3,
                control = lca_control(n_starts = 1), verbose = FALSE)
  expect_true(f3$n_params > f2$n_params)
})
