# Check if outcome is a Surv() expression

Tests whether an outcome specification string represents a survival
outcome by checking for the Surv() function pattern. Used to route model
fitting to Cox proportional hazards methods.

## Usage

``` r
is_surv_outcome(outcome)
```

## Arguments

- outcome:

  Character string of the outcome specification.

## Value

Logical `TRUE` if outcome starts with "Surv(", `FALSE` otherwise.
