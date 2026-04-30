# Apply locale decimal mark to a sprintf-formatted string

Replaces the period decimal mark in a `sprintf`-formatted string with
the locale-appropriate decimal mark, and fixes negative-zero artefacts.

## Usage

``` r
apply_decimal_mark(x, marks)
```

## Arguments

- x:

  Character string (already formatted with `sprintf`).

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.github.io/summata/reference/resolve_number_marks.md).

## Value

Character string with corrected decimal marks and no negative zeros.
