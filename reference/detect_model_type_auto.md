# Auto-detect model type based on outcome and random effects

Determines the appropriate model type based on outcome variable
characteristics and presence of random effects.

## Usage

``` r
detect_model_type_auto(data, outcome, has_random_effects, family = "binomial")
```

## Arguments

- data:

  Data.frame or data.table containing the outcome variable.

- outcome:

  Character string specifying the outcome variable or Surv() expression.

- has_random_effects:

  Logical indicating if random effects are specified.

- family:

  Character string for GLM family (default "binomial").

## Value

Character string indicating detected model type.
