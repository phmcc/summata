# Process survival probability quantiles (optimized)

Extracts survival time quantiles from survfit objects. Uses vectorized
operations for efficiency.

## Usage

``` r
process_survival_probs(
  survfit_objects,
  probs,
  groups,
  group_labels,
  time_digits,
  total,
  total_label,
  median_label,
  by,
  data,
  conf_level = 0.95,
  marks = NULL
)
```

## Arguments

- survfit_objects:

  List of survfit objects.

- probs:

  Numeric vector of probabilities.

- groups:

  Character vector of group names.

- group_labels:

  Character vector of group display labels.

- time_digits:

  Integer decimal places for time values.

- total:

  Logical or character controlling total column.

- total_label:

  Character label for total column.

- median_label:

  Character label for median row.

- by:

  Character name of stratifying variable.

- data:

  Data.table with the source data.

## Value

List with formatted and raw data.tables.
