# Extract metrics for mixed-effects Cox models (coxme)

Extracts quality metrics specific to mixed-effects Cox proportional
hazards models including concordance, pseudo-R-squared, and ICC.

## Usage

``` r
extract_coxme_metrics(model, raw_data, metrics)
```

## Arguments

- model:

  Fitted coxme object from coxme package.

- raw_data:

  Data.table with raw model information.

- metrics:

  Named list of initialized metrics to populate.

## Value

Updated metrics list with coxme-specific values.
