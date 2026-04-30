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
  marks = NULL
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

## Value

Formatted data.table.
