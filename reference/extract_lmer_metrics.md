# Extract metrics for linear mixed-effects models (lmer)

Extracts quality metrics specific to linear mixed-effects models
including R-squared measures, ICC, and global significance tests.

## Usage

``` r
extract_lmer_metrics(model, raw_data, metrics)
```

## Arguments

- model:

  Fitted lmerMod object from lme4.

- raw_data:

  Data.table with raw model information.

- metrics:

  Named list of initialized metrics to populate.

## Value

Updated metrics list with lmer-specific values.
