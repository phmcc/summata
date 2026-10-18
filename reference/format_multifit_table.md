# Format multifit results for publication

Internal helper that formats raw multivariate results into
publication-ready table format.

## Usage

``` r
format_multifit_table(
  data,
  columns,
  show_n = TRUE,
  show_events = TRUE,
  digits = 2,
  p_digits = 3,
  labels = NULL,
  predictor_label = NULL,
  include_predictor = TRUE,
  exponentiate = NULL,
  conf_level = 0.95,
  marks = NULL,
  interaction_vars = NULL,
  predictor_levels = NULL,
  interaction_levels = NULL
)
```

## Arguments

- data:

  Raw combined data.table from combine_multifit_results.

- columns:

  Character specifying column layout.

- show_n:

  Logical for sample size column.

- show_events:

  Logical for events column.

- digits:

  Integer decimal places for effects.

- p_digits:

  Integer decimal places for p-values.

- labels:

  Named vector of labels for outcomes and predictors.

- predictor_label:

  Label for predictor variable.

- include_predictor:

  Logical for including predictor column.

- exponentiate:

  Logical or `NULL` for coefficient handling.

- interaction_vars:

  Character vector of the variables of the interaction terms, by which
  their coefficients are labeled (see
  [`format_interaction_term()`](https://phmcc.codefloe.page/summata/reference/format_interaction_term.md)).

- predictor_levels:

  Character vector of the levels of the predictor, which are told apart
  from interaction terms when they contain a colon.

- interaction_levels:

  Named list of the levels of the variables in `interaction_vars`, which
  keeps a level that contains a colon in one component of an interaction
  term.

## Value

Formatted data.table.
