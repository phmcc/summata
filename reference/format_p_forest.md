# Format a p-value for forest plot annotations

Returns a formatted p-value string suitable for forest plot display,
using locale-aware decimal marks when `marks` is provided.

## Usage

``` r
format_p_forest(p, p_digits, marks = NULL)
```

## Arguments

- p:

  Numeric p-value.

- p_digits:

  Integer decimal places.

- marks:

  Optional list with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.codeberg.page/summata/reference/resolve_number_marks.md).

## Value

Character string.
