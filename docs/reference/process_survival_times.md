# Process survival at specified time points (optimized)

Extracts survival probabilities at specified time points from survfit
objects. Uses vectorized operations for efficiency.

## Usage

``` r
process_survival_times(
  survfit_objects,
  times,
  groups,
  group_labels,
  stats,
  type,
  digits,
  percent,
  total,
  total_label,
  time_label,
  time_unit,
  by,
  data,
  marks = NULL
)
```

## Arguments

- survfit_objects:

  List of survfit objects.

- times:

  Numeric vector of time points.

- groups:

  Character vector of group names.

- group_labels:

  Character vector of group display labels.

- stats:

  Character vector of statistics to include.

- type:

  Character string specifying probability type.

- digits:

  Integer decimal places for percentages.

- percent:

  Logical whether to display as percentages.

- total:

  Logical or character controlling total column.

- total_label:

  Character label for total column.

- time_label:

  Character template for time column headers.

- time_unit:

  Character time unit for column headers.

- by:

  Character name of stratifying variable.

- data:

  Data.table with the source data.

## Value

List with formatted and raw data.tables.
