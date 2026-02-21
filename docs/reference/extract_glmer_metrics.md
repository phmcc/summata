# Extract metrics for generalized linear mixed-effects models (glmer)

Extracts quality metrics specific to generalized linear mixed-effects
models including concordance, R-squared measures, ICC, and Brier score
for binomial.

## Usage

``` r
extract_glmer_metrics(model, raw_data, metrics)
```

## Arguments

- model:

  Fitted glmerMod object from lme4.

- raw_data:

  Data.table with raw model information.

- metrics:

  Named list of initialized metrics to populate.

## Value

Updated metrics list with glmer-specific values.
