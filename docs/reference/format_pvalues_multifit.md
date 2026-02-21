# Format p-values for multifit display

Converts numeric p-values to formatted character strings. Values below
the threshold (determined by digits parameter) display as "\< 0.001"
(for digits=3), "\< 0.0001" (for digits=4), *etc.* NA values display as
"-".

## Usage

``` r
format_pvalues_multifit(p, digits = 3, marks = NULL)
```

## Arguments

- p:

  Numeric vector of p-values.

- digits:

  Integer number of decimal places. Also determines the threshold for
  "less than" display: threshold = 10^(-digits). Default is 3.

- marks:

  Optional list with `big.mark` and `decimal.mark` as returned by
  [`resolve_number_marks`](https://phmcc.github.io/summata/reference/resolve_number_marks.md).

## Value

Character vector of formatted p-values.
