# Identify the coefficients of interactions

R joins the components of an interaction coefficient with colons
(*e.g.,* "treatmentDrug A:stageII"), but a factor level may contain a
colon as well (*e.g.,* the level "08:00" of a factor `shift`, whose
coefficient is "shift08:00"), as may the name of a variable. A name with
a colon is therefore taken for an interaction unless it is the
coefficient of a level in `xlevels` or the name of a main-effect term of
the model.

## Usage

``` r
is_interaction_term(terms, xlevels = NULL, main_terms = NULL)
```

## Arguments

- terms:

  Character vector of coefficient names, without the backticks of
  non-syntactic variable names.

- xlevels:

  Named list of factor levels, or `NULL`.

- main_terms:

  Character vector of the main-effect terms of the model, as returned by
  [`main_effect_labels()`](https://phmcc.codefloe.page/summata/reference/main_effect_labels.md),
  or `NULL`.

## Value

Logical vector, `TRUE` for the coefficients of interactions.
