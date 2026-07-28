# Identify the variables required by a fitted model

Returns the names of every variable referenced by a model formula,
including the response, stratification factors, offsets and
random-effect grouping factors. Terms are expanded where possible, so a
model specified with a dot on the right-hand side (*e.g.,* `y ~ .`)
returns the complete set of variables rather than the literal dot.

## Usage

``` r
get_model_variables(model, model_class)
```

## Arguments

- model:

  Fitted model object.

- model_class:

  Character string of the model class.

## Value

Character vector of variable names, or `character(0)` when the formula
cannot be recovered.
