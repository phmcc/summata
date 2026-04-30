# Format *p*-values for exported tables

Applies bold formatting to significant *p*-values in LaTeX tables using
the textbf command.

## Usage

``` r
format_pvalues_export_tex(df, p_threshold = 0.05)
```

## Arguments

- df:

  Data.table containing *p*-value columns.

- p_threshold:

  Numeric threshold for significance (default 0.05).

## Value

Data.table with significant *p*-values wrapped in LaTeX bold commands.
