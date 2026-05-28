## Pre-compute the slow fits used by the vignettes. Run with:
##   Rscript tools/build-vignette-artefacts.R
## Outputs go to inst/extdata/*.rds. These objects are loaded by the
## vignettes so the rendered docs show live numbers without paying
## the EM cost at vignette-build time (or at every devtools::check()).

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})
out_dir <- file.path("inst", "extdata")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Dataset 1: poLCA::election (categorical workflow vignette)
# ---------------------------------------------------------------------------
data("election", package = "poLCA")
cat_items <- c("MORALG","CARESG","KNOWG","LEADG","DISHONG","INTELG",
               "MORALB","CARESB","KNOWB","LEADB","DISHONB","INTELB")
keep <- stats::complete.cases(election[, cat_items])
elec <- election[keep, ]

cat("=== election: K = 2, 3, 4 (naive) ===\n")
elec_fits <- lapply(2:4, function(K) {
  fit_lca(elec, categorical = cat_items, n_classes = K,
          control = lca_control(n_starts = 5, seed = 110),
          verbose = FALSE)
})
names(elec_fits) <- paste0("K", 2:4)
saveRDS(elec_fits, file.path(out_dir, "election_naive_fits.rds"))

cat("=== election: K = 3 with BVR-guided direct effects ===\n")
elec_bvr <- auto_bvr(
  data = elec, categorical = cat_items,
  K_range = 3,                # fix K so the table reads cleanly
  max_direct_effects = 4L,
  bvr_threshold = 3.84,
  seed = 110, verbose = FALSE,
  n_starts = 5)
saveRDS(elec_bvr, file.path(out_dir, "election_bvr_fit.rds"))

cat("=== election: K = 3 with SLD (auto_sld) ===\n")
elec_sld <- auto_sld(
  data = elec, categorical = cat_items,
  n_classes = 3,
  max_rank_per_class = 3L,
  criterion = "BIC",
  seed = 110, verbose = FALSE)
saveRDS(elec_sld, file.path(out_dir, "election_sld_fit.rds"))

# ---------------------------------------------------------------------------
# Dataset 2: mlbench::PimaIndiansDiabetes2 (continuous + distal vignette)
# ---------------------------------------------------------------------------
data("PimaIndiansDiabetes2", package = "mlbench")
pima_vars <- c("glucose","pressure","mass","pedigree")
pima <- PimaIndiansDiabetes2[stats::complete.cases(
  PimaIndiansDiabetes2[, c(pima_vars,"age","diabetes")]), ]

cat("=== pima: K = 2, 3 (naive continuous) ===\n")
pima_naive <- lapply(2:3, function(K) {
  fit_lca(pima, continuous = pima_vars, n_classes = K,
          control = lca_control(n_starts = 10, seed = 110),
          verbose = FALSE)
})
names(pima_naive) <- paste0("K", 2:3)
saveRDS(pima_naive, file.path(out_dir, "pima_naive_fits.rds"))

cat("=== pima: K = 2 with concomitant age (character vec) ===\n")
pima_concom_chr <- fit_lca(
  pima, continuous = pima_vars, concomitant = "age",
  n_classes = 2,
  control = lca_control(n_starts = 10, seed = 110),
  verbose = FALSE)
saveRDS(pima_concom_chr, file.path(out_dir, "pima_concom_chr.rds"))

cat("=== pima: K = 2 with concomitant ~ age + I(age^2) (formula) ===\n")
pima_concom_fm <- fit_lca(
  pima, continuous = pima_vars, concomitant = ~ age + I(age^2),
  n_classes = 2,
  control = lca_control(n_starts = 10, seed = 110),
  verbose = FALSE)
saveRDS(pima_concom_fm, file.path(out_dir, "pima_concom_fm.rds"))

cat("=== pima: penalized covariance ===\n")
pima_pen <- fit_lca(
  pima, continuous = pima_vars, concomitant = "age",
  n_classes = 2, dependence = "penalized",
  control = lca_control(n_starts = 10, seed = 110),
  verbose = FALSE)
saveRDS(pima_pen, file.path(out_dir, "pima_penalized.rds"))

cat("=== pima: distal model for diabetes ===\n")
pima_distal <- distal(pima_concom_chr, pima,
                      formula = diabetes ~ age,
                      family  = "binomial")
saveRDS(pima_distal, file.path(out_dir, "pima_distal.rds"))

cat("\nAll artefacts written to ", out_dir, "\n", sep = "")
