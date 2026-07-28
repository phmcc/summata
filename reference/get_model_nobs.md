# Number of observations used to fit a model

Provides a single accessor for the fitted sample size across the model
classes supported by summata. Survival models record this quantity in
class-specific components rather than through
[`stats::nobs()`](https://rdrr.io/r/stats/nobs.html).

## Usage

``` r
get_model_nobs(model, model_class)
```

## Arguments

- model:

  Fitted model object.

- model_class:

  Character string of the model class.

## Value

Numeric scalar giving the number of observations used in fitting, or
`NA_real_` when it cannot be determined.
