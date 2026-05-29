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
# Dataset 1: voter_perceptions (categorical workflow vignette)
# ---------------------------------------------------------------------------
data("voter_perceptions", package = "mixLCA")
cat_items <- names(voter_perceptions)

cat("=== voter_perceptions: K = 2, 3, 4 (naive) ===\n")
voter_fits <- lapply(2:4, function(K) {
  fit_lca(voter_perceptions, categorical = cat_items, n_classes = K,
          control = lca_control(n_starts = 5),
          verbose = FALSE)
})
names(voter_fits) <- paste0("K", 2:4)
saveRDS(voter_fits, file.path(out_dir, "voter_naive_fits.rds"))

cat("=== voter_perceptions: K = 3 with BVR-guided direct effects ===\n")
voter_bvr <- auto_bvr(
  data = voter_perceptions, categorical = cat_items,
  K_range = 3,
  max_direct_effects = 4L,
  bvr_threshold = 3.84,
  verbose = FALSE,
  n_starts = 5)
saveRDS(voter_bvr, file.path(out_dir, "voter_bvr_fit.rds"))

cat("=== voter_perceptions: K = 3 with SLD (auto_sld) ===\n")
voter_sld <- auto_sld(
  data = voter_perceptions, categorical = cat_items,
  n_classes = 3,
  max_rank_per_class = 3L,
  criterion = "BIC",
  verbose = FALSE)
saveRDS(voter_sld, file.path(out_dir, "voter_sld_fit.rds"))

# ---------------------------------------------------------------------------
# Dataset 2: health_screening (continuous + distal vignette)
# ---------------------------------------------------------------------------
data("health_screening", package = "mixLCA")
hs_vars <- c("marker_1", "marker_2", "marker_3", "marker_4")

cat("=== health_screening: K = 2, 3 (naive continuous) ===\n")
hs_naive <- lapply(2:3, function(K) {
  fit_lca(health_screening, continuous = hs_vars, n_classes = K,
          control = lca_control(n_starts = 10),
          verbose = FALSE)
})
names(hs_naive) <- paste0("K", 2:3)
saveRDS(hs_naive, file.path(out_dir, "hs_naive_fits.rds"))

cat("=== health_screening: K = 2 with concomitant age (character vec) ===\n")
hs_concom_chr <- fit_lca(
  health_screening, continuous = hs_vars, concomitant = "age",
  n_classes = 2,
  control = lca_control(n_starts = 10),
  verbose = FALSE)
saveRDS(hs_concom_chr, file.path(out_dir, "hs_concom_chr.rds"))

cat("=== health_screening: K = 2 with concomitant ~ age + I(age^2) (formula) ===\n")
hs_concom_fm <- fit_lca(
  health_screening, continuous = hs_vars, concomitant = ~ age + I(age^2),
  n_classes = 2,
  control = lca_control(n_starts = 10),
  verbose = FALSE)
saveRDS(hs_concom_fm, file.path(out_dir, "hs_concom_fm.rds"))

cat("=== health_screening: penalized covariance ===\n")
hs_pen <- fit_lca(
  health_screening, continuous = hs_vars, concomitant = "age",
  n_classes = 2, dependence = "penalized",
  control = lca_control(n_starts = 10),
  verbose = FALSE)
saveRDS(hs_pen, file.path(out_dir, "hs_penalized.rds"))

cat("=== health_screening: distal model for outcome ===\n")
hs_distal <- distal(hs_concom_chr, health_screening,
                    formula = outcome ~ age,
                    family  = "binomial")
saveRDS(hs_distal, file.path(out_dir, "hs_distal.rds"))

# ---------------------------------------------------------------------------
# Remove the legacy artefacts so they don't accumulate
# ---------------------------------------------------------------------------
legacy <- c("election_naive_fits.rds", "election_bvr_fit.rds",
            "election_sld_fit.rds",
            "pima_naive_fits.rds", "pima_concom_chr.rds",
            "pima_concom_fm.rds", "pima_penalized.rds",
            "pima_distal.rds")
for (f in file.path(out_dir, legacy)) {
  if (file.exists(f)) {
    file.remove(f)
    cat("removed legacy:", f, "\n")
  }
}

cat("\nAll artefacts written to ", out_dir, "\n", sep = "")
