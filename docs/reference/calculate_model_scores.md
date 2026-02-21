# Calculate Composite Mean Scores (CMS) for model comparison

Computes composite Score based on weighted combination of model quality
metrics. Weights vary by model type to reflect academic consensus on
important metrics for each model class.

## Usage

``` r
calculate_model_scores(comparison, model_type, scoring_weights = NULL)
```

## Arguments

- comparison:

  Data.table with model comparison metrics.

- model_type:

  Character string indicating model type.

- scoring_weights:

  Optional named list of custom weights. If `NULL`, uses default weights
  for the model type.

## Value

Data.table with CMS column added, sorted by score.
