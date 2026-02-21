# Format column headers with n counts (HTML)

Creates HTML-formatted column headers with sample size counts displayed
below the column name using line breaks.

## Usage

``` r
format_column_headers_with_n_html(col_names, n_row_data)
```

## Arguments

- col_names:

  Character vector of column names.

- n_row_data:

  Named list or data.table row with n values for each column.

## Value

Character vector with HTML-formatted headers including n counts.
