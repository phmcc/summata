# Format survival median with CI for display

Formats a survival median estimate with confidence interval using
locale-aware decimal marks and safe CI separators. Used by
`process_survival` in descriptive tables.

## Usage

``` r
format_survival_ci(median, lower, upper, fmt_str, marks)
```

## Arguments

- median:

  Numeric median survival time.

- lower:

  Numeric lower CI bound.

- upper:

  Numeric upper CI bound.

- fmt_str:

  Character string `sprintf` format specification.

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.codeberg.page/summata/reference/resolve_number_marks.md).

## Value

Character string with formatted "median (lower-upper)".
