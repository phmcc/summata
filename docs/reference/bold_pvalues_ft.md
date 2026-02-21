# Bold significant *p*-values in DOCX

Applies bold formatting to significant *p*-values in flextable objects
by detecting values below threshold or "\< 0.001" patterns.

## Usage

``` r
bold_pvalues_ft(ft, df, p_threshold = 0.05)
```

## Arguments

- ft:

  flextable object.

- df:

  The source data.table.

- p_threshold:

  Numeric *p*-value threshold for significance.

## Value

Flextable object with significant *p*-values bolded.
