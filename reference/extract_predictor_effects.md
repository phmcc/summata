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
  terms_to_extract = NULL,
  conf_method = NULL,
  predictor_levels = NULL
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

- predictor_levels:

  Character vector of the levels of the predictor, by which its
  coefficients are named (`"FALSE"` and `"TRUE"` for a logical
  predictor, none for a numeric one). When `NULL`, the levels recorded
  by the model are used.

## Value

data.table with predictor effect information.
