# Format an integer for forest plot annotations

Applies thousands separator to an integer value. Respects locale marks
when provided.

## Usage

``` r
format_count_forest(x, marks = NULL)
```

## Arguments

- x:

  Integer value.

- marks:

  Optional list with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.codeberg.page/summata/reference/resolve_number_marks.md).

## Value

Character string.
