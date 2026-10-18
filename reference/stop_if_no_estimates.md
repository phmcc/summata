# Stop when a model has no coefficient to report

A model none of whose coefficients other than the intercept could be
estimated, such as a model whose only predictor is constant, leaves
[`m2dt()`](https://phmcc.codefloe.page/summata/reference/m2dt.md) no row
to report once aliased coefficients are omitted. The table cannot be
formatted, so the cause is reported instead.

## Usage

``` r
stop_if_no_estimates(raw_data)
```

## Arguments

- raw_data:

  Data table returned by
  [`m2dt()`](https://phmcc.codefloe.page/summata/reference/m2dt.md)
  without the intercept.

## Value

`NULL`, invisibly, when `raw_data` has rows.
