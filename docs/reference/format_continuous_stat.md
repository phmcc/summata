# Format continuous statistic for display

Converts numeric summary statistics into formatted display strings
following standard conventions (mean ± SD, median \[IQR\], range,
*etc.*).

## Usage

``` r
format_continuous_stat(stats, stat_type, fmt_str, marks)
```

## Arguments

- stats:

  Named list of numeric statistics (mean, sd, median, q1, q3, min, max).

- stat_type:

  Character string: "mean_sd", "median_iqr", "median_range", or "range".

- fmt_str:

  Character string format specification for sprintf.

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.github.io/summata/reference/resolve_number_marks.md).

## Value

Character string with formatted statistic.
