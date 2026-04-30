# Identify variable groups before indentation

Detects variable group boundaries by finding rows where Variable column
is non-empty. Returns row indices for each group for zebra stripe
application.

## Usage

``` r
identify_variable_groups(df)
```

## Arguments

- df:

  Data.table with Variable column.

## Value

List of integer vectors, each containing row indices for one variable
group.
