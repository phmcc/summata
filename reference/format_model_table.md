# Format model results for publication-ready display

Transforms raw model coefficient data into a formatted table suitable
for publication. Handles effect measure formatting (OR, HR, RR,
Estimate), confidence intervals, *p*-values, sample sizes, and variable
labels. Supports interaction terms and mixed-effects models.

## Usage

``` r
format_model_table(
  data,
  effect_col = NULL,
  digits = 2,
  p_digits = 3,
  labels = NULL,
  show_n = TRUE,
  show_events = TRUE,
  reference_label = "reference",
  exponentiate = NULL,
  conf_level = 0.95,
  marks = NULL
)
```

## Arguments

- data:

  Data.table containing raw model results with coefficient columns.

- effect_col:

  Optional character string specifying the effect column name. If
  `NULL`, auto-detects from OR, HR, RR, or Estimate columns.

- digits:

  Integer number of decimal places for effect estimates.

- p_digits:

  Integer number of decimal places for *p*-values.

- labels:

  Optional named character vector mapping variable names to display
  labels. Supports automatic labeling of interaction terms.

- show_n:

  Logical whether to include sample size column.

- show_events:

  Logical whether to include events column (ignored for linear models).

- reference_label:

  Character string to display for reference categories.

- exponentiate:

  Optional logical to force exponentiated (TRUE) or raw (FALSE)
  coefficient display. If `NULL`, uses existing columns.

## Value

Formatted data.table with publication-ready columns.
