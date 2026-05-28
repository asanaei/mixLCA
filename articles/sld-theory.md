# Spectral Local Dependence: Theory and Worked Example

## Why a new method

The local-independence assumption underwrites every classical LCA: given
class $`k`$, the manifest indicators $`Y_1, \ldots, Y_J`$ are
statistically independent. In practice this assumption frequently fails.
Two established responses exist.

**Direct effects** (Hagenaars 1988; Vermunt 1999) add a conditional
probability table $`P(Y_j \mid Y_{j'}, C = k)`$ for one or a few item
pairs $`(j, j')`$. The remedy is targeted and cheap when local
dependence sits in a small number of identifiable pairs, but it scales
as $`O(J^2)`$ candidate pairs to search and produces dense parameter
tables for any pair admitted.

**Latent trait extensions** (Magidson and Vermunt’s mixed mixture
models, Bartolucci’s hidden-Markov LCA, factor-mixture models) add a
within-class continuous latent variable. They handle distributed
dependence but break the discrete-class story and introduce numerical
integration.

**Spectral Local Dependence (SLD)** sits between the two. It introduces
a rank-$`d_k`$*low-dimensional shift* to the class-conditional log
response surface. The shift is sufficient to capture broad correlated
residuals, but it leaves the class labels discrete and the M-step
closed-form (an eigendecomposition).

This vignette gives the math, the identification argument, the
implementation in `mixLCA`, and a worked example.

## The setup

Let $`J`$ be the number of categorical items, with item $`j`$ taking
values in $`\{1, \ldots, K_j\}`$. Stack the **one-hot encoding** of all
items into a single binary vector
``` math
\mathbf{Z}_i \in \{0,1\}^C, \qquad C := \sum_{j=1}^J K_j.
```
For each respondent $`i`$, exactly one entry of each item’s
$`K_j`$-block of $`\mathbf{Z}_i`$ equals 1. Let
$`\boldsymbol{\pi}^{(k)} \in [0,1]^C`$ be the class-$`k`$ marginal
response probabilities (the standard LCA parameters), one block per item
summing to 1. Under local independence, the row-wise expectation of
$`\mathbf{Z}`$ conditional on class $`k`$ is $`\boldsymbol{\pi}^{(k)}`$
and the conditional variance has a known multinomial form.

Define the **class-conditional residual**
``` math
\mathbf{R}^{(k)}_i := \mathbf{Z}_i - \boldsymbol{\pi}^{(k)}.
```
Its second moment under the class is the **conditional Burt matrix**
``` math
\boldsymbol{\Sigma}^{(k)} :=
  \mathbb{E}\!\left[\mathbf{R}^{(k)}\mathbf{R}^{(k)\top} \mid C = k\right].
```
Under local independence, the off-block-diagonal entries of
$`\boldsymbol{\Sigma}^{(k)}`$ are zero: items are uncorrelated given
class. Local dependence shows up as nonzero off-block-diagonal entries.

## The spectral idea

If we modelled $`\boldsymbol{\Sigma}^{(k)}`$ as a *full* matrix, we
would introduce $`\binom{C}{2}`$ free off-diagonal parameters per class,
far more than the data support. SLD instead assumes that the
off-diagonal covariance is *low-rank*: there exist $`d_k \ll C`$
orthonormal vectors
$`\mathbf{V}^{(k)} = [\mathbf{v}^{(k)}_1, \ldots, \mathbf{v}^{(k)}_{d_k}]
\in \mathbb{R}^{C \times d_k}`$ such that
``` math
\boldsymbol{\Sigma}^{(k)} \approx
  \mathbf{V}^{(k)} \boldsymbol{\Lambda}^{(k)} \mathbf{V}^{(k)\top}
```
where
$`\boldsymbol{\Lambda}^{(k)} = \mathrm{diag}(\lambda_1^{(k)}, \ldots,
\lambda_{d_k}^{(k)})`$ collects the top $`d_k`$ eigenvalues of
$`\boldsymbol{\Sigma}^{(k)}`$.

In words: each class has its own dominant directions of residual
covariance, and we model only those. The columns of $`\mathbf{V}^{(k)}`$
are *spectral loadings*: they tell you which item categories drift
together when conditional independence fails.

## The shifted log-likelihood

We do not parameterize $`\boldsymbol{\Sigma}^{(k)}`$ directly. Instead
we modify the class-$`k`$ log-density. Let
$`\boldsymbol{\eta}^{(k)}_i := \log \boldsymbol{\pi}^{(k)}`$ be the
local-independence log-probability vector. SLD adds a **hollow
projection** to that vector:
``` math
\boldsymbol{\eta}^{(k)}_i := \log \boldsymbol{\pi}^{(k)}
  + \mathbf{A}^{\star (k)} \mathbf{R}^{(k)}_i,
\qquad
\mathbf{A}^{\star (k)} := \mathbf{M} \odot \left(
  \mathbf{V}^{(k)} \mathbf{V}^{(k)\top}
\right),
```
where $`\mathbf{M}`$ is a $`C \times C`$ “hollow” mask (1
off-block-diagonal, 0 on-block-diagonal) that prevents the projection
from touching the *within-item* part of the parameters. The within-item
part is already parameterized by $`\boldsymbol{\pi}^{(k)}`$ and
double-counting it would break identification.

The log-density per item block is then a softmax of the shifted
$`\boldsymbol{\eta}^{(k)}_i`$ restricted to that item’s columns. The
composite log-likelihood is the sum across blocks.

Note: $`\mathbf{A}^{\star (k)}`$ is symmetric and rank-$`d_k`$ (before
masking; the mask preserves the rank up to numerical tolerance). The
identification of $`\mathbf{V}^{(k)}`$ is up to orthogonal rotation
within its column span, the same indeterminacy as in factor analysis or
correspondence analysis.

## Estimation

EM with one twist:

**E-step**. Compute posteriors $`\gamma_{ik} \propto \pi_k \cdot
\mathrm{softmax}_J(\boldsymbol{\eta}^{(k)}_i)`$, as usual.

**M-step**. For each class $`k`$:

1.  Update marginal response probabilities $`\boldsymbol{\pi}^{(k)}`$
    from the weighted posterior counts (standard LCA M-step).
2.  Form the weighted conditional residual covariance
    ``` math
    \widehat{\boldsymbol{\Sigma}}^{(k)} =
      \frac{1}{\sum_i \gamma_{ik}}
      \sum_i \gamma_{ik} \, \mathbf{R}^{(k)}_i \mathbf{R}^{(k)\top}_i,
    ```
    where
    $`\mathbf{R}^{(k)}_i = \mathbf{Z}_i - \boldsymbol{\pi}^{(k)}`$.
3.  Compute the eigendecomposition of
    $`\mathbf{M} \odot \widehat{\boldsymbol{\Sigma}}^{(k)}`$, keep the
    top $`d_k`$ eigenvectors, and form
    $`\mathbf{A}^{\star (k)} = \mathbf{M} \odot
    \mathbf{V}^{(k)} \mathbf{V}^{(k)\top}`$.
4.  Use this $`\mathbf{A}^{\star (k)}`$ to construct the shifted
    log-density for the next E-step.

The M-step is closed-form (an eigendecomposition), so the per-iteration
cost is $`O(K \cdot C^3)`$, dominated by eigendecomposition, but with
$`C`$ rarely exceeding a few dozen this is cheap.

EM monotonicity for this scheme rests on a generalized M-step argument:
the eigendecomposition selects the best rank-$`d_k`$ projection in
Frobenius norm, which is also the rank-$`d_k`$ projection that maximises
the Q-function locally. `mixLCA` adds a step-halving safeguard in
[R/04_em_engine.R](https://asanaei.github.io/mixLCA/R/04_em_engine.R)
that backs off if a candidate update fails to non-decrease the observed
log-likelihood; in practice the safeguard fires rarely.

## Choosing the ranks $`d_1, \ldots, d_K`$

Three strategies.

**Fixed rank.** Pass `spectral_rank = 2L` (scalar, same rank for all
classes) or `spectral_rank = c(1L, 2L, 0L)` (length-$`K`$ vector).

**Pooled basis.** Set `spectral_pool = TRUE` to estimate a single
spectral basis from the class-weighted average of conditional Burt
matrices, applied uniformly to every class. This sacrifices class-
specific dependence structure for a sharper estimate of the dominant
directions, useful when sample size per class is small.

**Adaptive search.**
[`auto_sld()`](https://asanaei.github.io/mixLCA/reference/auto_sld.md)
runs a greedy forward search: at each step it identifies the class with
the largest unmodelled eigenvalue, tries adding one rank there, and
accepts only if BIC (or AIC, aBIC, ICL) improves. The search stops when
no class can improve or when a rank cap is reached. The result is
class-specific ranks matched to the data’s actual residual structure.

## Identification

A few points worth flagging.

- The marginal probability vector $`\boldsymbol{\pi}^{(k)}`$ is the
  standard LCA parameter; identification of these probabilities relies
  on the same conditions as classical LCA (sufficient items, sufficient
  class separation).
- The rank-$`d_k`$ projection $`\mathbf{A}^{\star (k)}`$ is identified
  up to orthogonal rotation within its column span, exactly as in PCA,
  factor analysis, or multiple correspondence analysis.
- The hollow mask $`\mathbf{M}`$ is essential. Without it the spectral
  shift would absorb information already encoded in
  $`\boldsymbol{\pi}^{(k)}`$ (the within-item probabilities), producing
  a non-identified surplus.
- Because the M-step is composite (we maximise the sum of marginal item
  softmaxes, not the full multinomial joint), BIC computed in the usual
  way is *naive*: it under-penalises the composite likelihood. For rank
  selection in publication settings, prefer cross-validation over BIC.
  `mixLCA` flags this in
  [`summary()`](https://rdrr.io/r/base/summary.html) output when
  `spectral_rank > 0`.

## Worked example: when SLD wins

In
[`vignette("workflow")`](https://asanaei.github.io/mixLCA/articles/workflow.md)
we see that on the `voter_perceptions` dataset SLD beats both naive LCA
and BVR-guided direct effects, even though the bivariate residuals are
individually large. That is the typical case when local dependence is
broad. Below is a smaller synthetic example that makes the comparison
even sharper.

We generate $`N = 500`$ observations from two classes, each with 8
binary items. In class 1, items 1-4 are positively correlated via a
shared latent factor; in class 2 items 5-8 are positively correlated via
a different shared factor. Direct effects would need to enumerate
$`\binom{4}{2} \cdot 2 = 12`$ pairs to absorb this structure; SLD
captures it with two rank-1 projections.

``` r

sim_one_class <- function(n, p_base, latent_loading) {
  # latent ~ N(0,1); item j response is 1 with prob expit(logit(p_base[j]) + latent_loading[j] * latent)
  L <- rnorm(n)
  out <- sapply(seq_along(p_base), function(j) {
    eta <- stats::qlogis(p_base[j]) + latent_loading[j] * L
    rbinom(n, 1, plogis(eta))
  })
  as.data.frame(out)
}

withr::with_seed(110, {
  c1 <- sim_one_class(250,
    p_base = c(0.7,0.7,0.7,0.7, 0.3,0.3,0.3,0.3),
    latent_loading = c(2,2,2,2, 0,0,0,0))
  c2 <- sim_one_class(250,
    p_base = c(0.3,0.3,0.3,0.3, 0.7,0.7,0.7,0.7),
    latent_loading = c(0,0,0,0, 2,2,2,2))
  df_sim <- rbind(c1, c2)
})
df_sim[] <- lapply(df_sim, function(x) factor(x, levels = c(0,1), labels = c("no","yes")))
names(df_sim) <- paste0("Q", 1:8)
str(df_sim)
#> 'data.frame':    500 obs. of  8 variables:
#>  $ Q1: Factor w/ 2 levels "no","yes": 1 2 2 2 2 2 1 2 1 2 ...
#>  $ Q2: Factor w/ 2 levels "no","yes": 2 2 2 2 2 2 2 2 1 2 ...
#>  $ Q3: Factor w/ 2 levels "no","yes": 2 2 1 2 2 2 1 2 1 2 ...
#>  $ Q4: Factor w/ 2 levels "no","yes": 2 2 2 1 2 2 2 2 1 2 ...
#>  $ Q5: Factor w/ 2 levels "no","yes": 1 2 2 1 1 1 1 1 1 1 ...
#>  $ Q6: Factor w/ 2 levels "no","yes": 1 1 1 2 1 1 1 1 2 1 ...
#>  $ Q7: Factor w/ 2 levels "no","yes": 2 1 1 1 2 1 2 2 1 2 ...
#>  $ Q8: Factor w/ 2 levels "no","yes": 1 1 1 1 1 2 1 1 2 1 ...
```

``` r

# Naive K=2
fit_naive_sim <- fit_lca(df_sim,
                         categorical = names(df_sim),
                         n_classes = 2,
                         control = lca_control(n_starts = 5, seed = 110),
                         verbose = FALSE)

# SLD K=2 with rank 1 per class
fit_sld_sim <- fit_lca(df_sim,
                       categorical = names(df_sim),
                       n_classes = 2,
                       spectral_rank = c(1L, 1L),
                       control = lca_control(n_starts = 5, seed = 110),
                       verbose = FALSE)

data.frame(
  Model    = c("Naive", "SLD (rank 1 each)"),
  log_lik  = round(c(fit_naive_sim$log_lik, fit_sld_sim$log_lik), 2),
  n_params = c(fit_naive_sim$n_params, fit_sld_sim$n_params),
  BIC      = round(c(fit_indices(fit_naive_sim)$BIC,
                     fit_indices(fit_sld_sim)$BIC), 2)
)
#>               Model  log_lik n_params     BIC
#> 1             Naive -2622.30       17 5350.25
#> 2 SLD (rank 1 each) -2497.49       31 5187.63
```

The SLD model captures the within-class correlation with one rank per
class and gains BIC despite the extra parameters (naive BIC
$`\approx 5350`$ vs. SLD BIC $`\approx 5188`$, a savings of about 160
points). Inspecting the loadings shows the rank-1 direction in class 1
as a bipolar contrast: items Q1-Q4 (the items whose latent factor was
active in class 1) load with one sign, items Q5-Q8 load with the
opposite sign. The single dimension thereby encodes the covariance
pattern induced by the shared latent factor:

``` r

get_loadings(fit_sld_sim)[
  get_loadings(fit_sld_sim)$class == "Class 1" &
    get_loadings(fit_sld_sim)$dimension == 1, ]
#>      class dimension item category     loading
#> 1  Class 1         1   Q1       no  0.21629277
#> 2  Class 1         1   Q1      yes -0.21629277
#> 3  Class 1         1   Q2       no  0.33012727
#> 4  Class 1         1   Q2      yes -0.33012727
#> 5  Class 1         1   Q3       no  0.24362194
#> 6  Class 1         1   Q3      yes -0.24362194
#> 7  Class 1         1   Q4       no  0.24315168
#> 8  Class 1         1   Q4      yes -0.24315168
#> 9  Class 1         1   Q5       no -0.02299066
#> 10 Class 1         1   Q5      yes  0.02299066
#> 11 Class 1         1   Q6       no -0.28537661
#> 12 Class 1         1   Q6      yes  0.28537661
#> 13 Class 1         1   Q7       no -0.30480933
#> 14 Class 1         1   Q7      yes  0.30480933
#> 15 Class 1         1   Q8       no -0.22557023
#> 16 Class 1         1   Q8      yes  0.22557023
```

## When SLD does **not** win

Two failure modes are worth knowing about.

**Pairwise dependence is concentrated in a few item pairs.** Then direct
effects beat SLD because they spend parameters only on the pairs that
need them; SLD has to allocate non-zero coefficients to every pair that
the latent direction projects onto.

**Sample size per class is small** relative to $`C`$. The conditional
Burt matrix is a $`C \times C`$ covariance estimate; estimating its
top-$`d_k`$ eigendecomposition requires enough class-weighted
observations to stabilise the eigenvectors. As a rough rule:
$`N_k \gtrsim 10 \cdot C / K`$.

In either case
[`auto_sld()`](https://asanaei.github.io/mixLCA/reference/auto_sld.md)
will tell you by returning `spectral_rank = 0` for the affected classes:
the BIC criterion will reject the spectral expansion.

## Connection to multiple correspondence analysis

Readers from the categorical-data tradition will notice the kinship with
multiple correspondence analysis (MCA): MCA computes eigendecompositions
of a Burt matrix to embed categorical items in a Euclidean space. SLD is
the *conditional* analogue, applied class-by-class within a mixture
model.

The connection is also computational: MCA’s leading dimensions recover
the dominant residual covariance structure of the *pooled* sample; SLD’s
per-class eigendecomposition recovers it within each latent class.

## References

- Bartholomew, D. J., Steele, F., Galbraith, J., & Moustaki, I. (2008).
  *Analysis of Multivariate Social Science Data*. Chapman & Hall.
- Hagenaars, J. A. (1988). Latent structure models with direct effects
  between indicators. *Sociological Methods & Research*, 16(3), 379-405.
- Magidson, J., & Vermunt, J. K. (2001). Latent class factor and cluster
  models, bi-plots, and related graphical displays. *Sociological
  Methodology*, 31, 223-264.
- Vermunt, J. K. (1999). A general nonparametric approach to the
  analysis of ordinal categorical data. *Sociological Methodology*, 29,
  187-223.

The methodological paper for SLD itself is in preparation; this vignette
will be updated with a citation once it is available.

## Session info

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] future_1.70.0 mixLCA_1.0.1 
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6        future.apply_1.20.2 jsonlite_2.0.0     
#>  [4] dplyr_1.2.1         compiler_4.6.0      tidyselect_1.2.1   
#>  [7] Rcpp_1.1.1-1.1      parallel_4.6.0      jquerylib_0.1.4    
#> [10] globals_0.19.1      systemfonts_1.3.2   scales_1.4.0       
#> [13] textshaping_1.0.5   yaml_2.3.12         fastmap_1.2.0      
#> [16] ggplot2_4.0.3       R6_2.6.1            generics_0.1.4     
#> [19] knitr_1.51          tibble_3.3.1        desc_1.4.3         
#> [22] bslib_0.11.0        pillar_1.11.1       RColorBrewer_1.1-3 
#> [25] rlang_1.2.0         cachem_1.1.0        xfun_0.57          
#> [28] fs_2.1.0            sass_0.4.10         S7_0.2.2           
#> [31] cli_3.6.6           pkgdown_2.2.0       withr_3.0.2        
#> [34] magrittr_2.0.5      digest_0.6.39       grid_4.6.0         
#> [37] lifecycle_1.0.5     vctrs_0.7.3         evaluate_1.0.5     
#> [40] glue_1.8.1          listenv_0.10.1      farver_2.1.2       
#> [43] codetools_0.2-20    ragg_1.5.2          parallelly_1.47.0  
#> [46] rmarkdown_2.31      tools_4.6.0         pkgconfig_2.0.3    
#> [49] htmltools_0.5.9
```
