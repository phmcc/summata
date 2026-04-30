# Format combined fullfit output from formatted tables

Merges univariable and multivariable results into a single
publication-ready table with side-by-side display. Uses vectorized merge
instead of per-variable loops.

## Usage

``` r
format_fullfit_combined(
  uni_formatted,
  multi_formatted,
  uni_raw,
  multi_raw,
  predictors,
  columns,
  metrics,
  show_n,
  show_events,
  labels,
  exponentiate = NULL,
  conf_level = 0.95
)
```

## Arguments

- uni_formatted:

  Formatted data.table from univariable screening.

- multi_formatted:

  Formatted data.table from multivariable model.

- uni_raw:

  Raw data.table with univariable coefficients.

- multi_raw:

  Raw data.table with multivariable coefficients.

- predictors:

  Character vector of all predictor names.

- columns:

  Character string specifying columns to show ("both", "uni", "multi").

- metrics:

  Character vector specifying metrics to show ("effect", "p").

- show_n:

  Logical whether to include sample size column.

- show_events:

  Logical whether to include events column.

- labels:

  Optional named character vector of variable labels.

- exponentiate:

  Optional logical for coefficient exponentiation.

- conf_level:

  Numeric confidence level for CI label.

## Value

Combined data.table with aligned univariable and multivariable results.
