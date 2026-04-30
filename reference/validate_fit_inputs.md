# Complete input validation for fit functions

Master validation function called by fit(), uniscreen(), fullfit().
Performs comprehensive checks on data structure, variable existence,
numeric parameter ranges, and model-outcome consistency. Returns
validated parameters with auto-corrections applied when appropriate.

## Usage

``` r
validate_fit_inputs(
  data,
  outcome,
  predictors,
  model_type,
  family = NULL,
  conf_level = 0.95,
  digits = 2,
  p_digits = 3,
  p_threshold = NULL,
  auto_correct_model = TRUE
)
```

## Arguments

- data:

  Data frame or data.table containing all variables.

- outcome:

  Character string outcome specification (may include Surv()).

- predictors:

  Character vector of predictor variable names.

- model_type:

  Character string model type to validate.

- family:

  GLM family object, function, or string if applicable.

- conf_level:

  Numeric confidence level (must be between 0 and 1).

- digits:

  Integer number of decimal places for effect estimates.

- p_digits:

  Integer number of decimal places for *p*-values.

- p_threshold:

  Numeric *p*-value threshold for significance highlighting.

- auto_correct_model:

  Logical whether to auto-correct model type mismatches.

## Value

List with validated model_type, family, auto_corrected flag.
