# Fix negative zero in formatted strings

Corrects floating-point rounding artifacts that produce "-0.00" or
similar negative zero strings. Works on character vectors, replacing
patterns like "-0.00", "-0.000", *etc.* with their positive equivalents,
even when embedded within larger strings (*e.g.,* "(-0.00, 1.23)"
becomes "(0.00, 1.23)").

## Usage

``` r
fix_negative_zero(x, marks = NULL)
```

## Arguments

- x:

  Character vector of formatted numbers.

- marks:

  Optional list with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.github.io/summata/reference/resolve_number_marks.md).
  When `NULL`, uses the default US period decimal.

## Value

Character vector with negative zeros corrected.

## Details

When `marks` is supplied, also replaces the period decimal mark with the
locale-appropriate decimal mark.
