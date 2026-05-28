# Extract next unmodeled eigenvalue per class

For each class with rank below the ceiling, retrieves the next
eigenvalue from the conditional Burt matrix spectrum. If the model
already stores the full spectrum (as it does after at least one EM cycle
with SLD active), the value is read directly; otherwise, the weighted
residual Burt matrix is computed from the posteriors.

## Usage

``` r
get_unmodeled_eigenvalues(model, data, ranks, max_rank_per_class)
```
