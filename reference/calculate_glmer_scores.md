# Calculate scores for generalized linear mixed-effects models

Computes component scores for glmer models based on concordance,
marginal R-squared, and ICC metrics.

## Usage

``` r
calculate_glmer_scores(comparison, weights, scores, n_models)
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
