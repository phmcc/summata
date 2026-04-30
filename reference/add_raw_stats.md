# Add raw statistics to row

Appends raw numeric statistics to a data.table row for downstream
processing. Used to preserve underlying values alongside formatted
display strings.

## Usage

``` r
add_raw_stats(row, col, stats, stat_type)
```

## Arguments

- row:

  Data.table row to modify.

- col:

  Character string column name for the statistics.

- stats:

  Named list of numeric statistics (mean, sd, median, *etc.*).

- stat_type:

  Character string indicating which statistic is displayed.

## Value

Modified row (by reference).
