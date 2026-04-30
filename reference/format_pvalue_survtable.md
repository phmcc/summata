# Format *p*-value for survtable

Provides *p*-value formatting to the survtable result.

## Usage

``` r
format_pvalue_survtable(p, digits, marks = NULL)
```

## Arguments

- p:

  Numeric *p*-value.

- digits:

  Integer decimal places.

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.codeberg.page/summata/reference/resolve_number_marks.md).

## Value

Character-formatted *p*-value.
