# Classes of the fixed-effect variables of a mixed model

The counterpart, for lme4 models, of the `"dataClasses"` attribute of
the terms of other models: the class of each variable of the
fixed-effect terms, as the model frame records it. A variable is an
expression of the formula, such as `log(age)`, rather than a name within
one, so that the variables are named as the coefficients of their terms
begin.

## Usage

``` r
fixed_effect_classes(model)
```

## Arguments

- model:

  Fitted lme4 model.

## Value

Named character vector of the classes of the fixed-effect variables,
without the response.
