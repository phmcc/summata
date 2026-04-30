# Combine multivariate results

Internal helper to combine unadjusted and/or adjusted results from
multiple outcomes into a single data.table.

## Usage

``` r
combine_multifit_results(all_results, columns)
```

## Arguments

- all_results:

  List of results from fit_one_outcome calls.

- columns:

  Character specifying which columns to include.

## Value

Combined data.table.
