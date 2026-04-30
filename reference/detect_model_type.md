# Detect if model is univariable or multivariable

Determines whether a model contains one predictor (univariable) or
multiple predictors (multivariable) by analyzing coefficient names and
factor structure. Handles interactions and random effects appropriately.

## Usage

``` r
detect_model_type(model)
```

## Arguments

- model:

  Fitted model object.

## Value

Character string: "Univariable" or "Multivariable".
