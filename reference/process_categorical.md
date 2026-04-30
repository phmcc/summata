# Process categorical variable

Calculates frequency and percentage statistics for categorical
variables, with optional grouping and chi-square/Fisher's exact testing.
Handles factor levels, missing values, and custom labeling.

## Usage

``` r
process_categorical(
  data,
  var,
  var_label,
  group_var,
  stats,
  na_include,
  na_label,
  test,
  test_type,
  total,
  total_label,
  na_percent,
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

- na_include:

  Logical whether to include missing values as a category.

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

- na_percent:

  Logical whether to include NA in percentage denominators.

- ...:

  Additional arguments passed to test functions.

## Value

List with 'formatted' and 'raw' data.table components.
