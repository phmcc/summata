# Format numeric value with fixed decimal places

Formats a numeric value to a specified number of decimal places,
removing leading/trailing whitespace and fixing negative zero display
(*e.g.,* "-0.00" becomes "0.00"). When `marks` is supplied, applies
locale-appropriate decimal mark substitution.

## Usage

``` r
format_number(x, digits, marks = NULL)
```

## Arguments

- x:

  Numeric value to format.

- digits:

  Integer number of decimal places.

- marks:

  Optional list with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.codeberg.page/summata/reference/resolve_number_marks.md).

## Value

Character string with formatted value.
