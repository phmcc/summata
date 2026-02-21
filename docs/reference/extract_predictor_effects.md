# Extract predictor effects from a fitted model

Internal helper function that extracts only the predictor variable's
coefficient(s) from a fitted model, ignoring intercept and covariates.
Supports standard models (glm, lm, coxph) and mixed effects models
(glmer, lmer, coxme).

## Usage

``` r
extract_predictor_effects(
  model,
  predictor,
  outcome,
  conf_level = 0.95,
  adjusted = FALSE,
  terms_to_extract = NULL
)
```

## Arguments

- model:

  Fitted model object.

- predictor:

  Character string of the predictor variable name.

- outcome:

  Character string of the outcome variable name.

- conf_level:

  Numeric confidence level.

- adjusted:

  Logical indicating if this is an adjusted model.

- terms_to_extract:

  Character vector of terms to extract (predictor and optionally
  interaction terms involving the predictor).

## Value

data.table with predictor effect information.
