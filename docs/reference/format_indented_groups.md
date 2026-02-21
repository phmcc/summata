# Apply formatting to indented groups

Transforms tables with Variable/Group columns into indented format where
group levels appear as indented rows under variable names. Handles both
regression and descriptive tables with appropriate *p*-value placement.

## Usage

``` r
format_indented_groups(df, indent_string = "    ")
```

## Arguments

- df:

  Data.table with Variable and Group columns.

- indent_string:

  Character string to use for indentation.

## Value

Data.table with Group column removed and levels indented under
Variables.
