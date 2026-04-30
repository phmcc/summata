# Process variable wrapper

Routes variable processing to appropriate handler based on variable type
(continuous, categorical, or survival). Returns both formatted display
strings and raw numeric values.

## Usage

``` r
process_variable(
  data,
  var,
  group_var = NULL,
  stats_continuous,
  stats_categorical,
  digits,
  conf_level = 0.95,
  na_include,
  na_label,
  test,
  test_continuous,
  test_categorical,
  total,
  total_label,
  labels,
  na_percent,
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

- group_var:

  Optional character string naming the grouping variable.

- stats_continuous:

  Character vector of statistics for continuous variables.

- stats_categorical:

  Character vector of statistics for categorical variables.

- digits:

  Integer number of decimal places for continuous statistics.

- conf_level:

  Numeric confidence level for survival confidence intervals.

- na_include:

  Logical whether to include missing values as a category.

- na_label:

  Character string label for missing values.

- test:

  Logical whether to perform statistical tests.

- test_continuous:

  Character string specifying test type for continuous variables.

- test_categorical:

  Character string specifying test type for categorical variables.

- total:

  Logical or character controlling total column display.

- total_label:

  Character string label for total column.

- labels:

  Named character vector of variable labels.

- na_percent:

  Logical whether to include NA in percentage denominators.

- p_per_stat:

  Logical whether to show separate *p*-values per statistic for
  continuous variables. Default `FALSE` for better performance.

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.github.io/summata/reference/resolve_number_marks.md).

- ...:

  Additional arguments passed to test functions.

## Value

List with 'formatted' and 'raw' data.table components.
