# Condense table rows for more compact display

Collapses multi-row variables into single rows for compact tables.
Continuous variables show only the first statistic row, binary
categorical variables show only the non-reference category, and survival
variables show only the median row.

## Usage

``` r
condense_table_rows(df, indent_groups = TRUE)
```

## Arguments

- df:

  Data.table with Variable and Group columns.

- indent_groups:

  Logical whether indentation will be applied (affects processing).

## Value

Data.table with condensed rows.
