# Format *p*-values for descriptive tables

Converts numeric p-values to formatted strings with appropriate
precision. Handles very small p-values with threshold notation (*e.g.,*
"\< 0.001").

## Usage

``` r
format_pvalues_desctable(result, p_digits, marks)
```

## Arguments

- result:

  Data.table with 'p_value' column to format.

- p_digits:

  Integer number of decimal places for p-values.

- marks:

  List with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.github.io/summata/reference/resolve_number_marks.md).

## Value

Modified data.table with 'p_value' column (formatted strings).
