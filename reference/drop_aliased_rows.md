# Omit the rows of aliased coefficients from a forest plot

Removes the rows that display aliased coefficients, which have no
estimate and would otherwise be drawn as reference categories, and the
rows of any factor none of whose levels remains estimated, which would
otherwise be shown as a reference category compared with nothing. The
regression tables omit both. A factor whose first level has no
observations is estimated against an aliased level instead (see
[`effective_reference_level()`](https://phmcc.codefloe.page/summata/reference/effective_reference_level.md)).
That level is kept, first among the rows of the factor, as its reference
category, and the empty first level is removed, as in the regression
tables.

## Usage

``` r
drop_aliased_rows(terms_df, aliased, coef_terms, model = NULL)
```

## Arguments

- terms_df:

  Data table with one row per displayed term, with columns `var`,
  `level`, `N` and `pos` and, for interaction rows, `.int_term`.

- aliased:

  Character vector of aliased coefficient names, as returned by
  [`aliased_coefficients()`](https://phmcc.codefloe.page/summata/reference/aliased_coefficients.md).

- coef_terms:

  Character vector of the names of the model's coefficients.

- model:

  Optional fitted model, from which the reference of a factor with
  several aliased levels is found (see
  [`factor_reference()`](https://phmcc.codefloe.page/summata/reference/factor_reference.md)).

## Value

`terms_df` without the omitted rows.
