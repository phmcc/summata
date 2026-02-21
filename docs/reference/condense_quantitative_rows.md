# Condense quantitative variable rows only

Collapses multi-row continuous and survival variables into single rows
while preserving all categorical variable rows (including binary). Only
applies to descriptive tables from desctable().

## Usage

``` r
condense_quantitative_rows(df, indent_groups = TRUE)
```

## Arguments

- df:

  Data.table or data frame

- indent_groups:

  Logical. Whether to apply indentation formatting.

## Value

A data.table with condensed continuous/survival rows
