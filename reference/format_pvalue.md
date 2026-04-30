# Format a *p*-value with locale-aware decimal mark

Converts a numeric *p*-value to a display string with the correct
decimal separator and threshold notation (*e.g.,* "\< 0.001" or "\<
0,001").

## Usage

``` r
format_pvalue(p, digits, marks)
```

## Arguments

- p:

  Numeric *p*-value.

- digits:

  Integer number of decimal places.

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.github.io/summata/reference/resolve_number_marks.md).

## Value

Character string with the formatted *p*-value.
