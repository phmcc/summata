# Identify the event indicator among a model's response variables

Determines which response variable carries the event count, and whether
the model has one at all. Linear and Gaussian models do not; the
families that count events are those for which
[`m2dt()`](https://phmcc.codeberg.page/summata/reference/m2dt.md)
reports an events figure, so that a denominator is only offered where a
numerator exists.

## Usage

``` r
get_event_variable_for_counts(outcome_vars, model_type = NULL, family = NULL)
```

## Arguments

- outcome_vars:

  Character vector of response variable names.

- model_type:

  Character string of the model type or class.

- family:

  Model family, where applicable. Accepted as a name, a family object,
  or a generator function, since callers resolve the family at different
  points:
  [`fit()`](https://phmcc.codeberg.page/summata/reference/fit.md) passes
  the name it was given, while
  [`uniscreen()`](https://phmcc.codeberg.page/summata/reference/uniscreen.md)
  may already have resolved `"Gamma"` to `Gamma(link = "log")`.

## Value

Character name of the event indicator, or `NULL` when the model carries
no event count.
