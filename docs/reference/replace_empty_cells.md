# Replace empty cells with "-"

Converts empty strings and NA values to "-" for consistent display in
exported tables. Preserves Variable column values.

## Usage

``` r
replace_empty_cells(df)
```

## Arguments

- df:

  Data.frame or data.table to process.

## Value

Data.table with empty cells replaced by "-".
