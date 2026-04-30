# Vectorized quantile cell formatting

Formats survival quantile cells for multiple rows at once.

## Usage

``` r
format_quantile_cells(est, lower, upper, fmt_str, marks = NULL)
```

## Arguments

- est:

  Numeric vector of estimates.

- lower:

  Numeric vector of lower CI bounds.

- upper:

  Numeric vector of upper CI bounds.

- fmt_str:

  Format string for numeric values.

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.codeberg.page/summata/reference/resolve_number_marks.md).

## Value

Character vector of formatted cells.
