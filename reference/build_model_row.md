# Build comparison row for successfully fitted model

Creates a comparison table row with extracted metrics for a successfully
fitted model.

## Usage

``` r
build_model_row(
  model_name,
  n_predictors,
  converged,
  metrics,
  model_type,
  marks = NULL
)
```

## Arguments

- model_name:

  Character string name of the model.

- n_predictors:

  Integer number of predictors in the model.

- converged:

  Character string convergence status.

- metrics:

  Named list of extracted model metrics.

- model_type:

  Character string indicating model type.

## Value

Data.table with single row of formatted metrics.
