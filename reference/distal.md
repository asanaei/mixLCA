# Distal Outcome Estimation via BCH Weighting

Estimates class-specific regression models for one or more distal
outcomes, using BCH inverse-classification-error weights. The
measurement model posteriors are fixed before this step (the "cut"): no
gradient from distal outcomes reaches the class definitions.

## Usage

``` r
distal(model, data, formula, family = "gaussian")
```

## Arguments

- model:

  Fitted `mixLCA` object.

- data:

  Data frame containing both the original variables and the distal
  outcome.

- formula:

  A formula for the distal model, e.g.  
  `outcome ~ predictor1 + predictor2`.  
  A right-hand side of `~ 1` estimates unconditional class means.

- family:

  Character: `"gaussian"`, `"binomial"`, or `"poisson"`.

## Value

An object of class `mixDistal` containing class-specific model
summaries.
