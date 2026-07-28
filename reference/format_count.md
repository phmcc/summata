# Format integer counts with locale-aware thousands separator

Formats integer values with a thousands separator for display in tables.
Values below 1000 are returned as plain character strings, since a
grouping mark only applies from four digits.

## Usage

``` r
format_count(n, marks = NULL, na = NA_character_)
```

## Arguments

- n:

  Numeric vector of counts.

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.codeberg.page/summata/reference/resolve_number_marks.md).
  Resolved from the global option when not supplied, so that callers
  without a locale of their own still produce output consistent with the
  rest of the package.

- na:

  Character string substituted for missing values. Contexts differ in
  what reads correctly: tables conventionally show a dash, plot
  annotations are better left blank, and the default preserves the
  missing value for callers doing their own substitution.

## Value

Character vector of formatted counts, the same length as `n`.

## Details

Vectorized over `n`. Each element is formatted independently, so counts
of differing magnitude are not padded to a common width.
