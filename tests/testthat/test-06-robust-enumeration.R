# ==============================================================================
# tests/testthat/test-06-robust-enumeration.R
# Tests for run_em_robust() and enumerate_lca().
# ==============================================================================

test_that("enumerate_lca() returns models for every K in k_range", {
  df  <- make_df_small()
  out <- enumerate_lca(
    data        = df,
    continuous  = c("x1","x2","x3"),
    k_range     = 2:3,
    n_starts    = 1,
    verbose     = FALSE
  )

  expect_type(out, "list")
  expect_true(all(c("models","comparison") %in% names(out)))
  expect_equal(length(out$models), 2L)
  expect_s3_class(out$models$K2, "mixLCA")
  expect_s3_class(out$models$K3, "mixLCA")
})

test_that("enumerate_lca() comparison table has K=2 and K=3 rows", {
  df  <- make_df_small()
  out <- enumerate_lca(
    data       = df,
    continuous = c("x1","x2","x3"),
    k_range    = 2:3,
    n_starts   = 1,
    verbose    = FALSE
  )
  tab <- out$comparison
  expect_equal(nrow(tab), 2L)
  expect_equal(sort(tab$K), c(2L, 3L))
})

test_that("enumerate_lca() BIC of K=2 is lower than K=3 for 2-class data", {
  # The data are generated from 2 classes, so K=2 should fit best by BIC
  df  <- make_df_small()
  out <- enumerate_lca(
    data       = df,
    continuous = c("x1","x2","x3"),
    k_range    = 2:3,
    n_starts   = 2,
    verbose    = FALSE
  )
  tab <- out$comparison
  bic2 <- tab$BIC[tab$K == 2]
  bic3 <- tab$BIC[tab$K == 3]
  expect_true(bic2 < bic3,
              label = "K=2 BIC lower than K=3 for 2-class data")
})

test_that("enumerate_lca() errors when validate_inputs fails", {
  df <- make_df_small()
  expect_error(
    enumerate_lca(df, continuous = "nonexistent", k_range = 2:3,
                  verbose = FALSE),
    regexp = "absent"
  )
})
