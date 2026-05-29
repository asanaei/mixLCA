# Predict Bridge for tidyclust (class probabilities)

Returns a data frame of posterior class probabilities with columns named
`.pred_Cluster_1`, `.pred_Cluster_2`, etc., as required by the tidyclust
predict interface.

## Usage

``` r
.lca_clust_predict_prob_mixLCA(object, new_data)
```

## Arguments

- object:

  A fitted `mixLCA` object.

- new_data:

  Data frame of new observations.

## Value

Data frame of posterior probabilities.
