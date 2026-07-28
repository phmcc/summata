# Identify the response variables of a fitted model

Returns the variables appearing on the left-hand side of a model
formula. A survival response contributes both the time and status
variables.

## Usage

``` r
get_outcome_variables(model)
```

## Arguments

- model:

  Fitted model object.

## Value

Character vector of response variable names, or `character(0)` when the
formula cannot be recovered.
