# Process continuous variable

Calculates descriptive statistics for continuous numeric variables, with
optional grouping and statistical testing. Supports multiple summary
statistics (mean \\\pm\\ SD, median \[IQR\], range) and various
hypothesis tests.

## Usage

``` r
process_continuous(
  data,
  var,
  var_label,
  group_var,
  stats,
  digits,
  na_include,
  na_label,
  test,
  test_type,
  total,
  total_label,
  p_per_stat = FALSE,
  marks = NULL,
  ...
)
```

## Arguments

- data:

  Data.table containing the variable.

- var:

  Character string naming the variable to process.

- var_label:

  Character string label for display.

- group_var:

  Optional character string naming the grouping variable.

- stats:

  Character vector of statistics to calculate.

- digits:

  Integer number of decimal places.

- na_include:

  Logical whether to include missing values.

- na_label:

  Character string label for missing values.

- test:

  Logical whether to perform statistical tests.

- test_type:

  Character string specifying test type.

- total:

  Logical or character controlling total column display.

- total_label:

  Character string label for total column.

- p_per_stat:

  Logical. If TRUE, calculate separate *p*-values for each statistic
  type (*e.g.,* t-test for means, Wilcoxon for medians). If `FALSE`
  (default), calculate a single *p*-value based on the first statistic
  type.

- ...:

  Additional arguments passed to test functions.

## Value

List with 'formatted' and 'raw' data.table components.
