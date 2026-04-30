# Resolve the CI or range separator

Determines the appropriate separator character between two numeric
bounds (*e.g.,* CI lower-upper, range min-max) based on whether either
bound is negative and the current locale's decimal mark. This avoids
ambiguous output like `"(-5--3)"` or `"(1,2-3,4)"` with European commas.

## Usage

``` r
resolve_separator(lower, upper, marks)
```

## Arguments

- lower:

  Numeric value of the lower bound (or minimum).

- upper:

  Numeric value of the upper bound (or maximum).

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.codeberg.page/summata/reference/resolve_number_marks.md).

## Value

Character string separator.

## Details

Rules:

1.  If either bound is negative, use `" to "` to avoid double-hyphen
    ambiguity (*e.g.,* `"-5 to -3"` not `"-5--3"`).

2.  If the decimal mark is a comma (EU locale), use `"\u2013"` (en-dash)
    to avoid confusion between decimal commas and separating commas
    (*e.g.,* `"1,2\u20133,4"` not `"1,2-3,4"`).

3.  Otherwise, use a plain hyphen `"-"`.
