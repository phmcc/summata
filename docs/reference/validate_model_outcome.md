# Validate model type matches outcome specification

Ensures consistency between the specified model type, outcome variable
type, and GLM family (if applicable). Detects common mismatches like
using survival outcomes with non-survival models or binary outcomes with
linear models. Can auto-correct fixable issues or raise informative
errors.

## Usage

``` r
validate_model_outcome(
  outcome,
  model_type,
  family = NULL,
  data = NULL,
  auto_correct = TRUE
)
```

## Arguments

- outcome:

  Character string outcome specification (may include Surv()).

- model_type:

  Character string specified model type.

- family:

  GLM family object, function, or string if applicable.

- data:

  Data frame or data.table for outcome type detection.

- auto_correct:

  Logical whether to auto-correct fixable mismatches.

## Value

List with model_type, family, messages, auto_corrected flag.

## Details

Checks for mismatches and auto-corrects or errors as appropriate.
