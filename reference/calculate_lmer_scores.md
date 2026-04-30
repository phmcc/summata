# Calculate scores for linear mixed-effects models

Computes component scores for lmer models based on marginal R-squared,
conditional R-squared, and ICC metrics.

## Usage

``` r
calculate_lmer_scores(comparison, weights, scores, n_models)
```

## Arguments

- comparison:

  Data.table with model comparison metrics.

- weights:

  Named list of scoring weights.

- scores:

  List of initialized score vectors.

- n_models:

  Integer number of models being compared.

## Value

Updated scores list with calculated values and total.
