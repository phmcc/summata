# Extract comprehensive model metrics based on academic consensus

Extracts quality control metrics from fitted models for comparison.
Supports GLM, Cox, linear, and mixed-effects models.

## Usage

``` r
extract_model_metrics(model, raw_data, model_type)
```

## Arguments

- model:

  Fitted model object.

- raw_data:

  Data.table with raw model information.

- model_type:

  Character string indicating model type.

## Value

Named list of metrics.
