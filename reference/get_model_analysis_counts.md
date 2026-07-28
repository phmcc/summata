# Count the observations available to a fitted model

Convenience wrapper around
[`get_analysis_counts()`](https://phmcc.codeberg.page/summata/reference/get_analysis_counts.md)
that derives the response and predictor variables from the model itself.

## Usage

``` r
get_model_analysis_counts(model, model_class, data)
```

## Arguments

- model:

  Fitted model object.

- model_class:

  Character string of the model class.

- data:

  Data frame or data.table supplied to the model call.

## Value

See
[`get_analysis_counts()`](https://phmcc.codeberg.page/summata/reference/get_analysis_counts.md).
