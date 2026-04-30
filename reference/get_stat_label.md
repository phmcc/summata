# Get display label for statistic type

Converts internal statistic type codes to formatted display labels for
table column headers.

## Usage

``` r
get_stat_label(stat_type)
```

## Arguments

- stat_type:

  Character string: "mean_sd", "median_iqr", "median_range", "range",
  "n_miss", or custom type.

## Value

Character string with formatted label (*e.g.,* "Mean \\\pm\\ SD").
