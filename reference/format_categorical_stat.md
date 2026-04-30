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

  Integer total count for percentage calculation.

- stat_type:

  Character string: `"n"`, `"n_percent"`, or `"percent"`.

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.codeberg.page/summata/reference/resolve_number_marks.md).

## Value

Character string with the formatted statistic.
