# Extract event variable from survival model

Parses the Surv() expression in survival model formulas to extract the
event/status variable name. Works with coxph, clogit, and coxme models.

## Usage

``` r
get_event_variable(model, model_class)
```

## Arguments

- model:

  Fitted survival model object.

- model_class:

  Character string of the model class.

## Value

Character string naming the event variable, or `NULL` if not found.
