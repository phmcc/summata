# Process a single survival outcome

Parses a Surv() expression, fits the survival models, performs the group
comparison test where requested, and assembles the formatted and raw
tables for one outcome. Called once per outcome by survtable().

## Usage

``` r
process_single_outcome(
  data,
  outcome,
  outcome_label,
  by,
  times,
  probs,
  stats,
  type,
  conf_level,
  conf_type,
  digits,
  time_digits,
  percent,
  test,
  test_type,
  total,
  total_label,
  time_unit,
  time_label,
  median_label,
  labels,
  na_rm,
  marks = NULL,
  ...
)
```

## Arguments

- data:

  Data.table with the source data.

- outcome:

  Character string giving a Surv() expression.

- outcome_label:

  Character label for the outcome, or NULL for a single-outcome table.

- by:

  Character name of stratifying variable.

- times:

  Numeric vector of time points.

- probs:

  Numeric vector of probabilities.

- stats:

  Character vector of statistics to include.

- type:

  Character string specifying probability type.

- conf_level:

  Numeric confidence level for interval estimates.

- conf_type:

  Character string specifying confidence interval type.

- digits:

  Integer decimal places for percentages.

- time_digits:

  Integer decimal places for time values.

- percent:

  Logical whether to display as percentages.

- test:

  Logical whether to perform a group comparison test.

- test_type:

  Character string specifying test type.

- total:

  Logical or character controlling total column.

- total_label:

  Character label for total column.

- time_unit:

  Character time unit for column headers.

- time_label:

  Character template for time column headers.

- median_label:

  Character label for median row.

- labels:

  Named character vector of display labels for group levels.

- na_rm:

  Logical whether to drop rows with missing time, status, or stratifying
  values.

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.codeberg.page/summata/reference/resolve_number_marks.md).

- ...:

  Additional arguments passed to
  [`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html).

## Value

List with the formatted table, the raw table, the survfit objects, and
the test result.
