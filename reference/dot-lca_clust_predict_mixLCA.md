# Predict Bridge for tidyclust (cluster assignments)

Returns a factor of cluster labels. tidyclust's `predict.cluster_fit`
wraps this into the standard `.pred_cluster` tibble column.

## Usage

``` r
.lca_clust_predict_mixLCA(object, new_data)
```

## Arguments

- object:

  A fitted `mixLCA` object.

- new_data:

  Data frame of new observations.

## Value

Factor of cluster assignments.
