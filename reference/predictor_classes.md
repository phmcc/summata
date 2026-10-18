# Read the classes of the predictors of a terms object

Returns the classes of the variables of the model other than the
response, leaving out expressions such as `I(age > 65)`, which are not
variables of the data and whose levels cannot be counted there.

## Usage

``` r
predictor_classes(model_terms)
```

## Arguments

- model_terms:

  Terms object of a fitted model, or `NULL`.

## Value

Named character vector of the `"dataClasses"` of those variables, or
`NULL`.
