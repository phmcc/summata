# Get readable model type name

Converts model class names to human-readable descriptions. For GLMs,
uses the family to provide specific names (*e.g.,* "Logistic",
"Poisson").

## Usage

``` r
get_model_type_name(model)
```

## Arguments

- model:

  Fitted model object.

## Value

Character string with readable model type name.
