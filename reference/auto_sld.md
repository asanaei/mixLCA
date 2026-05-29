# Automated Class-Specific Spectral Rank Selection

Evaluates the necessity of Spectral Local Dependence (SLD) iteratively.
Begins with a base model (all ranks zero or a user-supplied fitted
model) and conducts a forward-stepwise greedy search. At each step the
unmodeled eigenvalue spectrum is inspected per class, and the class
exhibiting the largest residual eigenvalue is targeted for a rank
increment. If the resulting model improves the information criterion,
the increment is accepted; otherwise the search terminates.

## Usage

``` r
auto_sld(
  data,
  continuous = NULL,
  categorical = NULL,
  concomitant = NULL,
  n_classes = 2L,
  max_rank_per_class = 3L,
  max_total_rank = 10L,
  criterion = c("BIC", "aBIC", "AIC", "ICL"),
  base_model = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  Data frame.

- continuous:

  Character vector or NULL.

- categorical:

  Character vector (at least 2 required).

- concomitant:

  Character vector or NULL.

- n_classes:

  Integer: number of latent classes.

- max_rank_per_class:

  Integer: max spectral rank for any class.

- max_total_rank:

  Integer: max sum of ranks across all classes.

- criterion:

  Character: information criterion to minimize (`"BIC"`, `"aBIC"`,
  `"AIC"`, or `"ICL"`).

- base_model:

  Optional pre-fitted `mixLCA` object. If NULL, an independence model (d
  = 0) is fitted automatically.

- verbose:

  Logical: print search trajectory.

- ...:

  Additional arguments forwarded to `fit_lca` (e.g., `max_iter`, `tol`,
  `dependence`).

## Value

A `mixLCA` object of the selected configuration, with an additional
element `auto_spectral_path` documenting the search history (data frame
with columns `step`, `class_incremented`, `ranks`, and the criterion
value).

## Details

Because each candidate is hot-started from the current model, the latent
class definitions remain anchored, preventing latent-class drift.
Classes that do not exhibit residual dependence retain `d = 0`, avoiding
unnecessary parameter proliferation.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(110)
data(voter_perceptions)
fit <- auto_sld(
  data        = voter_perceptions,
  categorical = names(voter_perceptions),
  n_classes   = 3,
  max_rank_per_class = 3L,
  criterion   = "BIC",
  verbose     = FALSE)
fit$specs$spectral_rank      # class-specific ranks selected
fit$auto_spectral_path        # accepted increments
fit_indices(fit)$BIC
} # }
```
