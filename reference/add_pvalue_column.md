# Add *p*-value column to result table

Adds formatted *p*-value column to the survtable result.

## Usage

``` r
add_pvalue_column(result, p_value, p_digits, marks = NULL)
```

## Arguments

- result:

  Data.table result.

- p_value:

  Numeric *p*-value.

- p_digits:

  Integer decimal places for *p*-value.

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.codeberg.page/summata/reference/resolve_number_marks.md).

## Value

Data.table with *p*-value column added.
