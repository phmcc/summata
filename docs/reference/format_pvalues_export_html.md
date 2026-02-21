# Format *p*-values for exported tables (HTML)

Applies bold formatting to significant *p*-values in HTML tables using
the b tag.

## Usage

``` r
format_pvalues_export_html(df, p_threshold = 0.05)
```

## Arguments

- df:

  Data.table containing *p*-value columns.

- p_threshold:

  Numeric threshold for significance (default 0.05).

## Value

Data.table with significant *p*-values wrapped in HTML bold tags.
