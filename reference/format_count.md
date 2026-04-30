# Format an integer count with locale-aware thousands separator

Formats integer values with a thousands separator for display in tables.
Values below 1000 are returned as plain character strings.

## Usage

``` r
format_count(n, marks)
```

## Arguments

- n:

  Integer count value.

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.github.io/summata/reference/resolve_number_marks.md).

## Value

Character string with the formatted count.
