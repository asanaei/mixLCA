# M-Step Update for Continuous Parameters

Updates class-specific mean and covariance using posterior weights.
Missing values are imputed through their conditional expectations under
the current parameters (the standard EM treatment of incomplete data).
Optionally applies L1 soft-thresholding to off-diagonal covariance
elements, producing data-driven sparsity.

## Usage

``` r
update_continuous(
  Y,
  weights,
  mu_old,
  Sigma_old,
  dependence = "full",
  penalty = 0
)
```

## Arguments

- Y:

  Numeric matrix (N x d), may contain `NA`.

- weights:

  Numeric vector of length N (posterior weights for this class).

- mu_old:

  Current mean vector (length d).

- Sigma_old:

  Current d x d covariance matrix.

- dependence:

  Character: `"none"`, `"full"`, or `"penalized"`.

- penalty:

  Numeric L1 penalty (used only when `dependence = "penalized"`).

## Value

List with elements `mean` (updated mean) and `covariance` (updated
covariance).
