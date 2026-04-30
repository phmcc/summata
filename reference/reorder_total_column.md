# Reorder columns to position total column

Rearranges data.table columns to place the total column in the specified
position (first, last, or default). Ensures proper ordering of Variable,
Group, total, group columns, and *p*-value.

## Usage

``` r
reorder_total_column(result, total, total_label)
```

## Arguments

- result:

  Data.table with columns to reorder.

- total:

  Logical or character: `TRUE`/"first" (total first), "last" (total
  last).

- total_label:

  Character string name of the total column.

## Value

Modified data.table with reordered columns.
