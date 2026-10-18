# Format a categorical statistic for display

Converts frequency counts into formatted display strings following
standard conventions (n, n (%), % only) with locale-aware decimal marks.

## Usage

``` r
format_categorical_stat(n, total, stat_type, marks)
```

## Arguments

- n:

  Integer count for the category.

- total:

  Integer total count for percentage calculation. When zero, as for a
  group without observations, the percentage is undefined: the count
  alone is returned, or an empty string for `"percent"`.

- stat_type:

  Character string: `"n"`, `"n_percent"`, or `"percent"`.

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.codefloe.page/summata/reference/resolve_number_marks.md).

## Value

Character string with the formatted statistic.
