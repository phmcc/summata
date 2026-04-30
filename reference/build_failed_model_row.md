# Build row for failed model

Creates a comparison table row with appropriate NA values for a model
that failed to fit.

## Usage

``` r
build_failed_model_row(model_name, n, n_predictors, model_type)
```

## Arguments

- model_name:

  Character string name of the model.

- n:

  Integer sample size.

- n_predictors:

  Integer number of predictors attempted.

- model_type:

  Character string indicating model type.

## Value

Data.table with single row of NA metrics and "Failed" convergence.
