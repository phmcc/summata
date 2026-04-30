# Vectorized survival cell formatting

Formats survival probability cells for multiple rows at once. Uses
locale-aware decimal marks and safe CI separators that avoid ambiguity
with negative values or decimal commas.

## Usage

``` r
format_survival_cells(
  est,
  lower,
  upper,
  n_risk,
  n_event,
  stats,
  fmt_est,
  fmt_ci_lower,
  fmt_ci_upper,
  percent,
  marks = NULL
)
```

## Arguments

- est:

  Numeric vector of estimates.

- lower:

  Numeric vector of lower CI bounds.

- upper:

  Numeric vector of upper CI bounds.

- n_risk:

  Integer vector of numbers at risk.

- n_event:

  Integer vector of event counts.

- stats:

  Character vector of statistics to include.

- fmt_est:

  Format string for estimate.

- fmt_ci_lower:

  Format string for lower CI bound.

- fmt_ci_upper:

  Format string for upper CI bound.

- percent:

  Logical whether percentages.

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.github.io/summata/reference/resolve_number_marks.md).

## Value

Character vector of formatted cells.
