# Validate outcome exists in data

Checks that the specified outcome variable (or survival variables within
Surv() expression) exists in the dataset. Raises informative error if
variables are missing. Handles both simple outcomes and Surv()
expressions.

## Usage

``` r
validate_outcome_exists(data, outcome)
```

## Arguments

- data:

  Data frame or data.table to check.

- outcome:

  Character string outcome specification (may include Surv()).

## Value

Invisible `TRUE` if validation passes, otherwise stops with error.
