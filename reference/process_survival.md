# Process survival variable

Calculates survival statistics including median survival times with
confidence intervals, with optional grouping and log-rank testing.
Parses Surv() expressions and uses survival package functions.

## Usage

``` r
process_survival(
  data,
  var,
  var_label,
  group_var,
  digits,
  conf_level = 0.95,
  na_include,
  na_label,
  test,
  total,
  total_label,
  marks = NULL,
  ...
)
```

## Arguments

- data:

  Data.table containing the survival variables.

- var:

  Character string with Surv() expression (*e.g.,* "Surv(time,
  status)").

- var_label:

  Character string label for display.

- group_var:

  Optional character string naming the grouping variable.

- digits:

  Integer number of decimal places.

- conf_level:

  Numeric confidence level for confidence intervals.

- na_include:

  Logical whether to include missing values.

- na_label:

  Character string label for missing values.

- test:

  Logical whether to perform log-rank test.

- total:

  Logical or character controlling total column display.

- total_label:

  Character string label for total column.

- ...:

  Additional arguments (currently unused).

## Value

List with 'formatted' and 'raw' data.table components.
