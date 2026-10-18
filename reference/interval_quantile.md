# Quantile for confidence intervals recomputed at another level

[`uniforest()`](https://phmcc.codefloe.page/summata/reference/uniforest.md)
and
[`multiforest()`](https://phmcc.codefloe.page/summata/reference/multiforest.md)
draw the intervals of their tables, and recompute them from the
coefficients and their standard errors when another confidence level is
requested. The tables hold the exact *t* intervals of linear models
(unless Wald intervals were requested) and record their residual degrees
of freedom; the other intervals are recomputed as Wald intervals.

## Usage

``` r
interval_quantile(conf_level, df_residual = NULL)
```

## Arguments

- conf_level:

  Numeric confidence level.

- df_residual:

  Numeric vector of residual degrees of freedom, one per row, missing
  where the table's intervals are not *t* intervals, or `NULL` when the
  table records none.

## Value

Numeric vector of quantiles, one per row, or a single normal quantile
when `df_residual` is `NULL`.
