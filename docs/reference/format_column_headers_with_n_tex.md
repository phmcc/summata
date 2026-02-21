# Format column headers with n counts (TeX)

Creates LaTeX-formatted column headers with sample size counts displayed
below the column name in a stacked format.

## Usage

``` r
format_column_headers_with_n_tex(col_names, n_row_data)
```

## Arguments

- col_names:

  Character vector of column names.

- n_row_data:

  Named list or data.table row with n values for each column.

## Value

Character vector with LaTeX-formatted headers including n counts.
