## ============================================================================
## data-raw/make-example-data.R
##
## Generates the two synthetic datasets shipped with mixLCA:
##   - data/voter_perceptions.rda
##   - data/health_screening.rda
##
## Both datasets are authored by the package author for didactic purposes.
## They are not derived from any real-world survey or clinical data, and are
## licensed under GPL (>= 3) along with the rest of the package.
##
## Run from the package root:
##   Rscript data-raw/make-example-data.R
## ============================================================================

stopifnot(requireNamespace("withr", quietly = TRUE))

## ----------------------------------------------------------------------------
## 1. voter_perceptions
##
## Synthetic 2000-subject perception study with 12 four-level Likert items
## rating two hypothetical political candidates ("A" and "B") on six
## attributes. The generative model is a three-class mixture:
##   Class 1 (40 %): strongly favours A, unfavourable toward B.
##   Class 2 (40 %): strongly favours B, unfavourable toward A.
##   Class 3 (20 %): moderate / ambivalent toward both candidates.
##
## Within-candidate item dependence is induced by a class-specific shared
## latent factor (loadings 1.2 on the six items of the favoured candidate),
## so BVR diagnostics light up on pairwise dependence and SLD recovers a
## rank-1 direction in the polarized classes.
##
## Substantive analog (not the source): rating-scale items in the American
## National Election Studies (ANES) where respondents rate presidential
## candidates on attribute dimensions.
## ----------------------------------------------------------------------------

LEVELS_LIKERT <- c("poor", "fair", "good", "excellent")
ATTRS         <- c("moral", "compet", "lead", "honest", "intel", "empath")

draw_likert <- function(n, base_logits, loading, latent) {
  ## Draw N four-level ordinal responses from a cumulative-logit model where
  ## a class-specific shared latent factor shifts the linear predictor.
  ##
  ## base_logits: length-3 vector of cumulative cut-points on the logit scale
  ##              (intercepts of P(Y <= 1), P(Y <= 2), P(Y <= 3)).
  ## loading:     scalar coefficient on the shared latent factor.
  ## latent:      length-N vector of latent factor scores (one per subject).
  ##
  ## Returns a length-N integer vector in 1:4.
  cum_eta <- outer(latent, base_logits, function(L, a) a - loading * L)
  cum_p   <- plogis(cum_eta)
  p1      <- cum_p[, 1]
  p2      <- cum_p[, 2] - p1
  p3      <- cum_p[, 3] - cum_p[, 2]
  p4      <- 1 - cum_p[, 3]
  P       <- cbind(p1, p2, p3, p4)
  P[P < 0] <- 0  # rare floating-point underflow
  P       <- P / rowSums(P)
  apply(P, 1, function(p) sample.int(4L, 1L, prob = p))
}

make_voter_perceptions <- function() {
  withr::with_seed(110, {
    N    <- 2000L
    cls  <- sample(1:3, N, replace = TRUE, prob = c(0.4, 0.4, 0.2))
    L    <- rnorm(N)  # shared latent factor (one per subject; loadings vary by class)

    cuts_favour <- c(-2.5, -1.0,  0.5)   # rate this candidate well
    cuts_neutral <- c(-1.0,  0.0,  1.0)  # neutral baseline
    cuts_against <- c( 0.5,  1.5,  2.5)  # rate this candidate poorly

    df <- data.frame(row.names = seq_len(N))

    for (a in ATTRS) {
      ## Items about candidate A
      cuts_A      <- ifelse(cls == 1, list(cuts_favour),
                     ifelse(cls == 2, list(cuts_against), list(cuts_neutral)))
      load_A      <- ifelse(cls == 1, 1.2, ifelse(cls == 2, 0.6, 0.4))

      ## Items about candidate B (mirror)
      cuts_B      <- ifelse(cls == 1, list(cuts_against),
                     ifelse(cls == 2, list(cuts_favour), list(cuts_neutral)))
      load_B      <- ifelse(cls == 1, 0.6, ifelse(cls == 2, 1.2, 0.4))

      a_codes <- integer(N)
      b_codes <- integer(N)
      for (k in 1:3) {
        idx <- which(cls == k)
        if (length(idx) == 0L) next
        a_codes[idx] <- draw_likert(
          length(idx),
          base_logits = switch(k, cuts_favour, cuts_against, cuts_neutral),
          loading     = switch(k, 1.2, 0.6, 0.4),
          latent      = L[idx])
        b_codes[idx] <- draw_likert(
          length(idx),
          base_logits = switch(k, cuts_against, cuts_favour, cuts_neutral),
          loading     = switch(k, 0.6, 1.2, 0.4),
          latent      = L[idx])
      }

      df[[paste0(a, "_A")]] <- factor(LEVELS_LIKERT[a_codes],
                                      levels = LEVELS_LIKERT, ordered = TRUE)
      df[[paste0(a, "_B")]] <- factor(LEVELS_LIKERT[b_codes],
                                      levels = LEVELS_LIKERT, ordered = TRUE)
    }
    df
  })
}

voter_perceptions <- make_voter_perceptions()
str(voter_perceptions)
cat("\nvoter_perceptions: nrow =", nrow(voter_perceptions),
    "  ncol =", ncol(voter_perceptions), "\n")
cat("Example cross-tab (moral_A x moral_B):\n")
print(table(voter_perceptions$moral_A, voter_perceptions$moral_B))


## ----------------------------------------------------------------------------
## 2. health_screening
##
## Synthetic 800-subject clinical-screening dataset. Two latent metabolic
## classes drive four continuous indicators; an antecedent age covariate
## shifts class membership; a downstream binary outcome depends on class
## (with mild age effect).
##
## Substantive analog (not the source): population-level metabolic screening
## studies that combine continuous biomarkers with follow-up disease status.
## We deliberately avoid using any real screening cohort because (a) license
## terms vary and (b) labelling participants by latent "metabolic class"
## reproduces a framing that we wish to leave to substantive researchers.
## ----------------------------------------------------------------------------

make_health_screening <- function() {
  withr::with_seed(110, {
    N <- 800L
    ## Age covariate (truncated normal, 25-70 range).
    age <- pmin(pmax(rnorm(N, mean = 45, sd = 12), 25), 70)

    ## Latent class shifts with age via multinomial logit:
    ##   logit(P(class = 2)) = -2.0 + 0.05 * (age - 45)
    p_class2 <- plogis(-2.0 + 0.05 * (age - 45))
    cls      <- ifelse(runif(N) < p_class2, 2L, 1L)

    ## Continuous indicators on the log scale (loosely log-normal in raw form).
    sigma <- 0.35
    df <- data.frame(
      marker_1 = exp(ifelse(cls == 1, log(80),  log(140)) + rnorm(N, 0, sigma)),
      marker_2 = exp(ifelse(cls == 1, log(70),  log(95))  + rnorm(N, 0, sigma)),
      marker_3 = exp(ifelse(cls == 1, log(22),  log(33))  + rnorm(N, 0, sigma * 0.7)),
      marker_4 = exp(ifelse(cls == 1, log(0.30), log(0.60)) + rnorm(N, 0, sigma)),
      age      = round(age, 1)
    )

    ## Distal outcome: positive screening probability depends on class with
    ## modest additional age effect.
    p_outcome <- plogis(ifelse(cls == 1, -2.0, 0.5) + 0.03 * (age - 45))
    df$outcome <- factor(ifelse(runif(N) < p_outcome, "yes", "no"),
                         levels = c("no", "yes"))

    df
  })
}

health_screening <- make_health_screening()
str(health_screening)
cat("\nhealth_screening: nrow =", nrow(health_screening),
    "  ncol =", ncol(health_screening), "\n")
cat("Outcome distribution:\n"); print(table(health_screening$outcome))


## ----------------------------------------------------------------------------
## Save to data/ with xz compression (smallest .rda format).
## ----------------------------------------------------------------------------
save(voter_perceptions, file = "data/voter_perceptions.rda",
     compress = "xz", compression_level = 9L)
save(health_screening, file = "data/health_screening.rda",
     compress = "xz", compression_level = 9L)

cat("\nSizes:\n")
print(file.info(c("data/voter_perceptions.rda", "data/health_screening.rda"))[, "size", drop = FALSE])
