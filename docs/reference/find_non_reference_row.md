# Find non-reference row for binary variable condensing

Identifies the non-reference row in a binary categorical variable by
checking for NA estimates (reference rows have NA). This is more robust
than assuming row position or matching specific strings like "Yes" or
"Positive".

## Usage

``` r
find_non_reference_row(var_rows, estimate_col = "estimate")
```

## Arguments

- var_rows:

  Data.table containing rows for a single variable.

- estimate_col:

  Character string naming the estimate column (*e.g.,* "estimate",
  "coef"). Default is "estimate".

## Value

Integer index of the non-reference row within var_rows, or `NULL` if
cannot be determined (*e.g.,* no NA estimates found, or multiple non-NA
rows).
