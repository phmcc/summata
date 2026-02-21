# Validate predictors exist in data

Checks that all specified predictor variables exist in the dataset.
Handles interaction terms (splits on ":"), mixed-effects random effects
(ignores "\|" syntax), and raises informative errors for missing
variables.

## Usage

``` r
validate_predictors_exist(data, predictors)
```

## Arguments

- data:

  Data frame or data.table to check.

- predictors:

  Character vector of predictor variable names.

## Value

Invisible `TRUE` if validation passes, otherwise stops with error.
