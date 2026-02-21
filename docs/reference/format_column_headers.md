# Apply formatting to column headers in exported tables (PDF/LaTeX)

Formats column headers for LaTeX output by escaping special characters,
italicizing 'n' and 'p', and optionally adding vertical spacing.

## Usage

``` r
format_column_headers(col_names, add_header_space = TRUE)
```

## Arguments

- col_names:

  Character vector of column names.

- add_header_space:

  Logical whether to add vertical padding.

## Value

Character vector with LaTeX-formatted column names.
