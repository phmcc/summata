# Count the observations available to an analysis

Summarizes how many of the supplied observations could enter a model,
and separates those excluded for a missing response from those excluded
for a missing covariate. The distinction matters: a missing response and
a missing covariate raise different questions about whether a
complete-case analysis is appropriate.

## Usage

``` r
get_analysis_counts(
  data,
  outcome_vars,
  predictor_vars = NULL,
  event_var = NULL
)
```

## Arguments

- data:

  Data frame or data.table supplied to the analysis.

- outcome_vars:

  Character vector of response variable names.

- predictor_vars:

  Character vector of predictor variable names, or a list of such
  vectors when several models are described.

- event_var:

  Character name of the event indicator, or `NULL` for models that carry
  no event count. Factors are coded as the modeling functions code them,
  with any level beyond the first counting as an event.

## Value

Named list with `n_supplied`, `n_analyzed`, `n_missing_outcome` and
`n_missing_predictor`, and, where an event indicator was supplied,
`events_supplied` and `events_analyzed`. `NULL` when the counts cannot
be determined.

## Details

Several models may be described at once by passing a list of predictor
sets, as when a univariable screen fits one model per predictor. The
analyzed count is then returned as one element per set.
