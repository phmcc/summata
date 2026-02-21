# Build comprehensive comparison table

Selects and renames columns based on model type following academic
consensus for reporting model fit statistics.

## Usage

``` r
build_comparison_table(comparison, model_type)
```

## Arguments

- comparison:

  Data.table with raw comparison metrics.

- model_type:

  Character string indicating model type.

## Value

Formatted data.table with appropriate columns for the model type.
