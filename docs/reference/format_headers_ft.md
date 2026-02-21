# Format headers for flextable

Applies formatting to flextable headers including italicizing 'n',
adding sample size counts from N row data, and bolding all headers.

## Usage

``` r
format_headers_ft(ft, has_n_row, n_row_data)
```

## Arguments

- ft:

  flextable object.

- has_n_row:

  Logical whether source data had an N row.

- n_row_data:

  Data from the N row for adding counts to headers.

## Value

Formatted flextable object.
