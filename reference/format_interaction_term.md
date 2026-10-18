# Format interaction term for display

Converts R's internal interaction term format (*e.g.,* "treatmentDrug
A:stageII") to a more readable format (*e.g.,* "Treatment (Drug A) ×
Stage (II)").

## Usage

``` r
format_interaction_term(term, labels = NULL, var_names = NULL, xlevels = NULL)
```

## Arguments

- term:

  Character string of the interaction term from model coefficients.

- labels:

  Optional named vector of labels for variable names.

- var_names:

  Optional character vector of the variables of the interaction. When
  given, each component is split into its variable and level by the
  longest variable name it begins with (see
  [`split_interaction_term()`](https://phmcc.codefloe.page/summata/reference/split_interaction_term.md)),
  which holds for any variable name; otherwise the variable name is
  taken to be the leading run of lowercase letters, digits, and
  underscores.

- xlevels:

  Optional named list of the levels of the variables in `var_names`,
  which keeps a level that contains a colon in one component.

## Value

Formatted interaction term string.
