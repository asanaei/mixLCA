# ==============================================================================
# tests/testthat/helper-mixLCA.R
# Shared fixtures loaded automatically by testthat before every test file.
# ==============================================================================

library(MASS)

# ---- Minimal two-class dataset (N=300) ----
make_df_small <- function(seed = 110) {
  set.seed(seed)
  N  <- 300
  cl <- ifelse(runif(N) < 0.5, 1L, 2L)
  data.frame(
    x1  = ifelse(cl == 1L, rnorm(N, 10, 2), rnorm(N, 4, 2)),
    x2  = ifelse(cl == 1L, rnorm(N, 10, 2), rnorm(N, 4, 2)),
    x3  = ifelse(cl == 1L, rnorm(N,  2, 1), rnorm(N, 8, 1)),
    age = rnorm(N, 40, 10),
    cat1 = ifelse(cl == 1L,
                  sample(c("A","B","C"), N, replace = TRUE,
                         prob = c(.6,.3,.1)),
                  sample(c("A","B","C"), N, replace = TRUE,
                         prob = c(.1,.3,.6))),
    outcome = ifelse(cl == 1L,
                     rnorm(N, 80, 10),
                     rnorm(N, 40, 10)),
    stringsAsFactors = FALSE
  )
}

# ---- Full dataset used in the demonstration (N=900) ----
make_df_full <- function(seed = 110) {
  set.seed(seed)
  N <- 900
  Age    <- rnorm(N, 40, 15)
  Status <- rnorm(N, 120, 25)
  X_true <- cbind(1, Age, Status)
  lp     <- X_true %*% c(-2.5, 0.06, -0.01)
  true_class <- ifelse(runif(N) < 1 / (1 + exp(-lp)), 2L, 1L)

  Metric_X <- Metric_Y <- Metric_Z <- numeric(N)
  Cat_Type <- character(N)
  for (i in seq_len(N)) {
    if (true_class[i] == 1L) {
      d <- MASS::mvrnorm(1, c(25,25,8),
                   matrix(c(6,4,0,4,6,0,0,0,3),3,3))
      Cat_Type[i] <- sample(c("A","B","C"), 1,
                             prob = c(.7,.2,.1))
    } else {
      d <- MASS::mvrnorm(1, c(12,12,18), diag(3,3))
      Cat_Type[i] <- sample(c("A","B","C"), 1,
                             prob = c(.1,.3,.6))
    }
    Metric_X[i] <- d[1]; Metric_Y[i] <- d[2]; Metric_Z[i] <- d[3]
  }
  Metric_Y[sample(N, 40)] <- NA
  data.frame(
    Metric_X, Metric_Y, Metric_Z, Cat_Type,
    Age, Status,
    Clinical_Outcome = ifelse(
      true_class == 1L,
      50 + 0.8*Age + rnorm(N,0,8),
      20 + 0.1*Age + rnorm(N,0,8))
  )
}
