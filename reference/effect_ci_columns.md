# Find the effect columns of a regression table

Regression tables name each effect column after its measure and the
confidence level of its intervals: "OR (95% CI)", "aHR (90% CI)", "Adj.
Coefficient (95% CI)", "Exp(Coef) (99% CI)". Descriptive and survival
tables have no such column; their interval columns, such as "Median (95%
CI)" in
[`survtable()`](https://phmcc.codefloe.page/summata/reference/survtable.md),
name a statistic rather than an effect measure.

## Usage

``` r
effect_ci_columns(col_names)
```

## Arguments

- col_names:

  Character vector of column names.

## Value

Character vector of the names that denote effect columns.
